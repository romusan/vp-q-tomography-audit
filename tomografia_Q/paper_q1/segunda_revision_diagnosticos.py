"""Event-block bootstrap and a controlled 2x2 mask/velocity experiment.

All outputs stay in paper_q1/segunda_revision. Archived campaigns are read only.
Bootstrap conditions on fixed archived predictions; it does not refit models.
Controlled fits reproduce the original unconstrained augmented least squares.
"""
import ast,json,math,shutil
from pathlib import Path
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.sparse.linalg import lsmr
R=Path(__file__).resolve().parent;ROOT=R.parent;OUT=R/'segunda_revision';OUT.mkdir(exist_ok=True)
for f in ('manuscrito.tex','manuscrito_gji.pdf','cifras.tex','cifras_revision.tex'):
 target=OUT/('antes_'+f)
 if not target.exists():shutil.copy2(R/f,target)
SEED=20261005;B=10000
C={'Reference':'sgc_final_q_extra_station_tomography','Expanded':'sgc_q_ampliado_300km','RayDomain':'sgc_q_mascara_rayos','Adopted':'sgc_q_mascara_rayos_damp400'}
frames={k:pd.read_csv(ROOT/'output_V2'/c/'data/sgc_1000_q_predictions.csv') for k,c in C.items()}
def erms(d,p):return float(100*(1-np.sqrt(np.sum((d-p)**2)/np.sum(d*d))))
stats={};blocks={}
for k,p in frames.items():
 terms=pd.read_csv(ROOT/'output_V2'/C[k]/'data/sgc_1000_q_station_terms.csv').set_index('station').station_term_s
 d=p.delta_tstar_s.to_numpy();pred=p.pred_final_s.to_numpy();u=(p.station_i.map(terms)-p.station_j.map(terms)).to_numpy();c=pred-u
 stats[k]={'events':int(p.event_id.nunique()),'pairs':len(p),'E_RMS':erms(d,pred),'structure_ratio':float(np.std(c)/np.std(d)),'station_ratio':float(np.std(u)/np.std(d)),'component_corr':float(np.corrcoef(c,u)[0,1]),'prediction_data_corr':float(np.corrcoef(pred,d)[0,1]),'data_rms_s':float(np.sqrt(np.mean(d*d))),'initial_residual_rms_s':float(np.sqrt(np.mean((d-p.pred_initial_s.to_numpy())**2))),'final_rms_s':float(np.sqrt(np.mean((d-pred)**2))),'residual_prediction_inner_product_mean':float(np.mean((d-pred)*pred))}
 blocks[k]=pd.DataFrame({'event':p.event_id,'sse':(d-pred)**2,'sd2':d*d}).groupby('event',sort=True).sum()
for k in ('RayDomain','Adopted'):
 assert frames[k][['event_id','station_i','station_j','delta_tstar_s']].equals(frames['Expanded'][['event_id','station_i','station_j','delta_tstar_s']])
draws={};rng=np.random.default_rng(SEED)
for labels in (['Reference'],['Expanded','RayDomain','Adopted']):
 n=len(blocks[labels[0]]);arrays=[blocks[k][['sse','sd2']].to_numpy() for k in labels]
 assert all(blocks[k].index.equals(blocks[labels[0]].index) for k in labels)
 samples=[[] for _ in labels]
 for start in range(0,B,250):
  weights=rng.multinomial(n,np.full(n,1/n),size=min(250,B-start))
  for arr,collect in zip(arrays,samples):
   summed=weights@arr;collect.append(100*(1-np.sqrt(summed[:,0]/summed[:,1])))
 for k,values in zip(labels,samples):
  draws[k]=np.concatenate(values);stats[k]['conditional_event_bootstrap_95CI']=np.percentile(draws[k],[2.5,97.5]).tolist()
contrasts={}
for hi,lo in [('RayDomain','Expanded'),('Adopted','RayDomain')]:
 v=draws[hi]-draws[lo];contrasts[f'{hi}_minus_{lo}']={'difference_percentage_points':stats[hi]['E_RMS']-stats[lo]['E_RMS'],'paired_event_bootstrap_95CI':np.percentile(v,[2.5,97.5]).tolist()}
np.savez_compressed(OUT/'bootstrap_draws.npz',**draws)
result={'bootstrap':{'replicates':B,'seed':SEED,'conditional_on_fixed_predictions':True,'sampling_unit':'event; all retained pairs within an event remain together','campaigns':stats,'paired_contrasts':contrasts}}
(OUT/'diagnosticos.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
print('Event bootstrap complete:',json.dumps(contrasts),flush=True)
# Extract the exact archived Laplacian implementation, without waveform imports.
tree=ast.parse((ROOT/'run_sgc_final_vp_q_preprocess.py').read_text(encoding='utf-8'))
selected=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in ('idx_to_kji','kji_to_idx','laplacian')]
ns={'np':np,'sp':sp};exec(compile(ast.Module(body=selected,type_ignores=[]),'<archived-laplacian>','exec'),ns)
base=ROOT/'output_V2'/C['Expanded'];newbase=ROOT/'output_V2'/C['RayDomain']
g=np.load(base/'models/sgc_1000_grid.npz');axes=[g[k] for k in ('xg','yg','zg')]
newg=np.load(newbase/'models/sgc_1000_grid.npz');assert all(np.array_equal(a,newg[k]) for a,k in zip(axes,('xg','yg','zg')))
velocities={key:np.load(folder/'models/sgc_1000_vp_final.npy') for key,folder in [('InheritedVp',base),('UpdatedVp',newbase)]}
masks={key:np.load(folder/'models/sgc_1000_q_active_mask.npy').astype(bool) for key,folder in [('InheritedMask',base),('RayMask',newbase)]}
p=frames['Expanded'];n=len(p);shape=velocities['InheritedVp'].shape;size=np.prod(shape)
def xy(lat,lon,z):return np.array([(lon+73)*111.32*np.cos(np.radians(7)),(lat-7)*111.32,z])
raycache={}
def ray(src,rec):
 dist=float(np.linalg.norm(rec-src));nstep=max(2,math.ceil(dist/2.5)+1);ds=dist/(nstep-1)
 xyz=src+np.linspace(0,1,nstep)[:,None]*(rec-src);inside=np.ones(nstep,bool)
 for i,ax in enumerate(axes):inside&=(xyz[:,i]>=ax.min())&(xyz[:,i]<=ax.max())
 xyz=xyz[inside];indices=[np.clip(np.rint((xyz[:,i]-ax[0])/(ax[1]-ax[0])).astype(int),0,len(ax)-1) for i,ax in enumerate(axes)]
 ix,iy,iz=indices;flat=(iz*len(axes[1])+iy)*len(axes[0])+ix
 ids,inverse=np.unique(flat,return_inverse=True)
 weights=[np.bincount(inverse,weights=ds/v[iz,iy,ix],minlength=len(ids)) for v in velocities.values()]
 return ids,weights
rowids=[];colids=[];vals=[[],[]]
for number,row in enumerate(p.itertuples(index=False)):
 parts=[]
 for side in ('i','j'):
  key=(row.event_id,getattr(row,'station_'+side))
  if key not in raycache:raycache[key]=ray(xy(row.event_lat,row.event_lon,row.event_depth_km),xy(getattr(row,'station_'+side+'_lat'),getattr(row,'station_'+side+'_lon'),0))
  parts.append(raycache[key])
 ids=np.union1d(parts[0][0],parts[1][0]);vv=np.zeros((2,len(ids)))
 for sign,(indices,weights) in zip((1,-1),parts):
  ix=np.searchsorted(ids,indices)
  for k,w in enumerate(weights):vv[k,ix]+=sign*w
 keep=np.any(vv!=0,axis=0);ids=ids[keep];vv=vv[:,keep]
 rowids.extend([number]*len(ids));colids.extend(ids.tolist())
 for k in range(2):vals[k].extend(vv[k].tolist())
 if (number+1)%10000==0:print('Kernel rows',number+1,flush=True)
matrices={key:sp.csr_matrix((values,(rowids,colids)),shape=(n,size)) for key,values in zip(velocities,vals)}
del rowids,colids,vals,raycache
stations=sorted(set(p.station_i)|set(p.station_j));positions={s:i for i,s in enumerate(stations)}
station_indices=np.column_stack([p.station_i.map(positions),p.station_j.map(positions)]).ravel()
T=sp.csr_matrix((np.tile([1.,-1.],n),(np.repeat(np.arange(n),2),station_indices)),shape=(n,len(stations)))
q0=np.broadcast_to(np.where(axes[2][:,None,None]<30,1/40,1/60),shape).ravel();d=p.delta_tstar_s.to_numpy()
experiments=[]
for vkey,G in matrices.items():
 initial=G@q0
 for mkey,mask in masks.items():
  active=np.flatnonzero(mask);Acell=G[:,active];L=ns['laplacian'](mask);nst=len(stations)
  A=sp.hstack([Acell,T],format='csr')
  aug=sp.vstack([A,sp.hstack([3*L,sp.csr_matrix((L.shape[0],nst))]),sp.hstack([.05*sp.eye(len(active)),sp.csr_matrix((len(active),nst))]),sp.hstack([sp.csr_matrix((nst,len(active))),.15*sp.eye(nst)])],format='csr')
  rhs=np.concatenate([d-initial,np.zeros(L.shape[0]+len(active)+nst)])
  sol=lsmr(aug,rhs,maxiter=1200);pred=initial+A@sol[0]
  row={'velocity':vkey,'mask':mkey,'observations':n,'active_nodes':len(active),'iterations':int(sol[2]),'istop':int(sol[1]),'E_RMS':erms(d,pred),'residual_rms_s':float(np.sqrt(np.mean((d-pred)**2)))}
  key=f'{vkey}_{mkey}'
  pd.DataFrame({'event_id':p.event_id,'delta_tstar_s':d,'prediction_linear_s':pred}).to_csv(OUT/(key+'_predictions.csv'),index=False)
  if (vkey,mkey) in [('InheritedVp','InheritedMask'),('UpdatedVp','RayMask')]:
   label='Expanded' if vkey=='InheritedVp' else 'RayDomain'
   row['archived_prediction_max_abs_difference_s']=float(np.max(np.abs(pred-frames[label].pred_final_s)))
   row['archived_residual_rms_difference_s']=row['residual_rms_s']-stats[label]['final_rms_s']
   # Floating-point summation order changes LSMR's stopping iteration.
   # Check pointwise predictions and aggregate misfit at explicit tolerances.
   assert row['archived_prediction_max_abs_difference_s']<1e-4 and abs(row['archived_residual_rms_difference_s'])<1e-7,row
  experiments.append(row);print('Controlled fit',json.dumps(row),flush=True)
result['controlled_experiment']={'station_damping':.15,'cell_smoothing':3.,'cell_damping':.05,'fixed_pairs':n,'fits':experiments,'interpretation_scope':'Unconstrained linear predictions only; not validation of clipped maps.'}
result['bounds']={}
for label,camp in [('Adopted',C['Adopted']),('HigherCeiling','sgc_q_techo1500')]:
 path=ROOT/'output_V2'/camp;q=np.load(path/'models/sgc_1000_q_pilot.npy');a=np.load(path/'models/sgc_1000_q_active_mask.npy').astype(bool)
 result['bounds'][label]={'active_nodes':int(a.sum()),'at_lower_bound':int((a&np.isclose(q,20,rtol=0,atol=1e-8)).sum()),'lower_bound_fraction_percent':float(100*(a&np.isclose(q,20,rtol=0,atol=1e-8)).sum()/a.sum()),'maximum_depth_active_km':float(max(z for z,layer in zip(axes[2],a) if layer.any()))}
(OUT/'diagnosticos.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
print('Second review diagnostics completed',flush=True)
