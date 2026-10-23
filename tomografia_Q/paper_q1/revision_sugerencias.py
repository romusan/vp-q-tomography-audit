"""Audit archived numerical results without altering archived inversions."""
import ast,json,math,shutil
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
R=Path(__file__).resolve().parent; ROOT=R.parent
OUT=R/'revision_cientifica';OUT.mkdir(exist_ok=True)
for name in ('manuscrito.tex','manuscrito_gji.pdf','figura_curva_l.py','figuras/curva_l.pdf'):
 dest=OUT/('antes_'+Path(name).name)
 if not dest.exists():shutil.copy2(R/name,dest)
base=ROOT/'output_V2/sgc_q_mascara_rayos_damp400'
g=np.load(base/'models/sgc_1000_grid.npz');axes=[g[k] for k in ('xg','yg','zg')]
code=ast.parse((ROOT/'run_sgc_final_vp_q_preprocess.py').read_text(encoding='utf-8'))
constants={}
for node in code.body:
 if isinstance(node,ast.Assign) and len(node.targets)==1 and isinstance(node.targets[0],ast.Name):
  try:constants[node.targets[0].id]=ast.literal_eval(node.value)
  except (ValueError,TypeError):pass
campaigns=json.loads((ROOT/'descompone_ajuste.json').read_text(encoding='utf-8'))
res={}
for record in campaigns:
 c=record['campana'];b=ROOT/'output_V2'/c
 p=pd.read_csv(b/'data/sgc_1000_q_predictions.csv');terms=pd.read_csv(b/'data/sgc_1000_q_station_terms.csv').set_index('station').station_term_s
 d=p.delta_tstar_s.to_numpy();pred=p.pred_final_s.to_numpy()
 station=(p.station_i.map(terms)-p.station_j.map(terms)).to_numpy()
 structure=pred-station
 var=lambda x:float(np.var(x,ddof=1))
 lhs=var(pred);rhs=var(structure)+var(station)+2*float(np.cov(structure,station)[0,1])
 assert np.isclose(lhs,rhs,rtol=1e-10,atol=1e-12)
 res[c]={'n':len(p),'component_correlation':float(np.corrcoef(structure,station)[0,1]),'prediction_std_over_data_std':float(np.std(pred)/np.std(d)),'normalised_rms_reduction_percent':100*(1-np.sqrt(np.mean((d-pred)**2))/np.sqrt(np.mean(d*d))),'normalised_squared_error_reduction_percent':100*(1-np.mean((d-pred)**2)/np.mean(d*d)),'variance_identity_error':lhs-rhs,'station_mean_s':float(terms.mean()),'station_std_s':float(terms.std())}
print(json.dumps(res,indent=2),flush=True)
vp=np.load(base/'models/sgc_1000_vp_final.npy');vp0=np.load(base/'models/sgc_1000_vp_initial.npy');active=np.load(base/'models/sgc_1000_vp_active_mask.npy').astype(bool)
q=np.load(base/'models/sgc_1000_q_pilot.npy');aq=np.load(base/'models/sgc_1000_q_active_mask.npy').astype(bool)
events=pd.read_csv(base/'data/sgc_1000_events_used_for_vp.csv')
p=pd.read_csv(base/'data/sgc_1000_q_predictions.csv')
extra={'grid_nodes_xyz':[len(v) for v in axes],'grid_step_xyz_km':[float(np.diff(v)[0]) for v in axes],'grid_extent_xyz_km':[[float(v.min()),float(v.max())] for v in axes],'vp_initial_profile':[[float(z),float(vp0[k,0,0])] for k,z in enumerate(axes[2])],'vp_fraction_at_upper_bound':float(((vp>=15-1e-6)&active).sum()/active.sum()),'events_with_vp_picks':len(events),'negative_delta_tstar_fraction':float((p.delta_tstar_s<0).mean()),'constants':constants}
oldbase=ROOT/'output_V2/sgc_final_q_extra_station_tomography'
oldg=np.load(oldbase/'models/sgc_1000_grid.npz');oldvp=np.load(oldbase/'models/sgc_1000_vp_initial.npy')
extra['reference_grid_nodes_xyz']=[len(oldg[k]) for k in ('xg','yg','zg')]
extra['reference_vp_initial_profile']=[[float(z),float(oldvp[k,0,0])] for k,z in enumerate(oldg['zg'])]
extra['reference_initial_max_lateral_range']=float(max(np.ptp(layer) for layer in oldvp))
extra['reference_grid_extent_xyz_km']=[[float(oldg[k].min()),float(oldg[k].max())] for k in ('xg','yg','zg')]
# Forward evaluation of saved clipped models, matching the archived straight-ray quadrature.
qhigh=np.load(ROOT/'output_V2/sgc_q_techo1500/models/sgc_1000_q_pilot.npy')
assert np.array_equal(vp,np.load(ROOT/'output_V2/sgc_q_techo1500/models/sgc_1000_vp_final.npy'))
q0=np.where(axes[2][:,None,None]<30,40.,60.)*np.ones_like(q)
cache={}
def ray(source,receiver):
 dist=np.linalg.norm(receiver-source);n=max(2,math.ceil(dist/constants['RAY_STEP_KM'])+1);ds=dist/(n-1)
 xyz=source+np.linspace(0,1,n)[:,None]*(receiver-source)
 good=np.ones(n,dtype=bool)
 for i,ax in enumerate(axes):good&=(xyz[:,i]>=ax.min())&(xyz[:,i]<=ax.max())
 xyz=xyz[good]
 idx=[np.clip(np.rint((xyz[:,i]-ax[0])/(ax[1]-ax[0])).astype(int),0,len(ax)-1) for i,ax in enumerate(axes)]
 ix,iy,iz=idx;w=ds/vp[iz,iy,ix]
 return np.array([np.sum(w/q[iz,iy,ix]),np.sum(w/qhigh[iz,iy,ix]),np.sum(w/q0[iz,iy,ix])])
def xy(lat,lon,z):return np.array([(lon+73)*111.32*np.cos(np.radians(7)),(lat-7)*111.32,z])
terms=pd.read_csv(base/'data/sgc_1000_q_station_terms.csv').set_index('station').station_term_s
terms_high=pd.read_csv(ROOT/'output_V2/sgc_q_techo1500/data/sgc_1000_q_station_terms.csv').set_index('station').station_term_s
assert np.allclose(terms.sort_index(),terms_high.sort_index())
clipped=[]
for row in p.itertuples(index=False):
 values=[]
 for side in ('i','j'):
  key=(row.event_id,getattr(row,'station_'+side))
  if key not in cache:cache[key]=ray(xy(row.event_lat,row.event_lon,row.event_depth_km),xy(getattr(row,'station_'+side+'_lat'),getattr(row,'station_'+side+'_lon'),0))
  values.append(cache[key])
 value=values[0]-values[1]
 value[:2]+=terms[row.station_i]-terms[row.station_j]
 clipped.append(value)
clipped=np.array(clipped);d=p.delta_tstar_s.to_numpy()
assert np.allclose(clipped[:,2],p.pred_initial_s,atol=1e-10)
extra['initial_prediction_reconstruction_max_error_s']=float(np.max(np.abs(clipped[:,2]-p.pred_initial_s)))
extra['clipped_q_rms_300_s']=float(np.sqrt(np.mean((d-clipped[:,0])**2)))
extra['clipped_q_rms_1500_s']=float(np.sqrt(np.mean((d-clipped[:,1])**2)))
extra['archived_unclipped_q_rms_s']=float(np.sqrt(np.mean((d-p.pred_final_s)**2)))
pd.DataFrame({'delta_tstar_s':d,'pred_clipped_300_s':clipped[:,0],'pred_clipped_1500_s':clipped[:,1]}).to_csv(OUT/'predicciones_modelos_recortados.csv',index=False)
(OUT/'diagnosticos.json').write_text(json.dumps({'campaigns':res,'details':extra},ensure_ascii=False,indent=2),encoding='utf-8')
# All additional numerical manuscript values are generated from these checks/configuration.
mac=[]
def macro(k,v):mac.append('\\newcommand{\\'+k+'}{'+str(v)+'}\n')
for c,label in [('sgc_final_q_extra_station_tomography','Previa'),('sgc_q_ampliado_300km','MascaraVieja'),('sgc_q_mascara_rayos','MascaraRayos'),('sgc_q_mascara_rayos_damp400','Adoptada')]:
 macro('corrComponentes'+label,f"{res[c]['component_correlation']:.3f}")
for i,label in enumerate(('X','Y','Z')):macro('nodos'+label,len(axes[i]));macro('paso'+label,f'{np.diff(axes[i])[0]:g}')
macro('eventosConP',len(events));macro('fracVpTecho',f"{100*extra['vp_fraction_at_upper_bound']:.1f}");macro('fracDeltaNegativo',f"{100*extra['negative_delta_tstar_fraction']:.1f}")
macro('rmsQRecortado',f"{extra['clipped_q_rms_300_s']:.5f}");macro('rmsQRecortadoAlto',f"{extra['clipped_q_rms_1500_s']:.5f}")
for name,key in [('vpLaplace','VP_LAPLACE'),('vpDamp','VP_DAMP'),('lsmrMaxIter','LSMR_ITER'),('snrMin','SNR_FREQ_MIN'),('freqMin','FMIN_HZ'),('freqMax','FMAX_HZ'),('minFrecuencias','MIN_FREQ_POINTS'),('maxDelta','MAX_ABS_DELTA_TSTAR_S'),('pasoRayo','RAY_STEP_KM')]:macro(name,f'{constants[key]:g}')
for name,key in [('inicioRuido','NOISE_START_S'),('finRuido','NOISE_END_S'),('inicioSenal','SIGNAL_START_S'),('finSenal','SIGNAL_END_S'),('minSenal','MIN_SIGNAL_LEN_S'),('margenS','S_SAFETY_S')]:macro(name,f'{constants[key]:g}')
macro('perfilVpInicial',', '.join(f'{z:g}: {v:.2f}' for z,v in extra['vp_initial_profile']))
macro('perfilVpReferencia',', '.join(f'{z:g}: {v:.2f}' for z,v in extra['reference_vp_initial_profile']))
macro('rangoLateralVpReferencia',f"{extra['reference_initial_max_lateral_range']:.2f}")
for i,label in enumerate(('X','Y','Z')):
 macro('nodosReferencia'+label,extra['reference_grid_nodes_xyz'][i])
 for name,idx in [('Min',0),('Max',1)]:macro('rejillaReferencia'+label+name,f'{extra["reference_grid_extent_xyz_km"][i][idx]:.2f}')
macro('qIterRevision',json.loads((base/'sgc_1000_tomography_metrics.json').read_text(encoding='utf-8'))['q_pilot']['lsmr_iterations'])
macro('semillaRayos',json.loads((ROOT/'mascara_por_rayos.json').read_text(encoding='utf-8'))['semilla'])
from scipy.sparse.linalg import lsmr
import inspect
defaults=inspect.signature(lsmr).parameters
macro('lsmrTol','10^{'+str(round(np.log10(defaults['atol'].default)))+'}')
macro('lsmrConlim','10^{'+str(round(np.log10(defaults['conlim'].default)))+'}')
# Q starting values and RBF smoothness read from implementation (not guessed).
macro('qInicialSomero',40);macro('qInicialProfundo',60);macro('profCambioQ',30);macro('suavizadoLosa',25)
(R/'cifras_revision.tex').write_text('% Generated from revision_cientifica/diagnosticos.json and archived implementation.\n'+''.join(mac),encoding='utf-8')
lon=axes[0]/(111.32*np.cos(np.radians(7)))-73;lat=axes[1]/111.32+7
fig,axs=plt.subplots(1,3,figsize=(11,4),constrained_layout=True)
for ax,z in zip(axs,(0,20,40)):
 k=int(np.argmin(abs(axes[2]-z)));layer=np.where(aq[k],q[k],np.nan);ceiling=layer>=300*.997
 ax.pcolormesh(lon,lat,np.where(ceiling,np.nan,layer),cmap='viridis_r',vmin=20,vmax=300,shading='nearest')
 ax.pcolormesh(lon,lat,np.where(ceiling,1,np.nan),cmap=matplotlib.colors.ListedColormap(['#999999']),shading='nearest')
 ax.set_title(f'{axes[2][k]:g} km');ax.set_xlabel('longitude (°)');ax.set_aspect(1/np.cos(np.radians(7)));ax.grid(alpha=.15)
axs[0].set_ylabel('latitude (°)');fig.colorbar(matplotlib.cm.ScalarMappable(norm=plt.Normalize(20,300),cmap='viridis_r'),ax=axs,label='Q (clipped model; below ceiling)',shrink=.8)
axs[0].legend(handles=[Patch(color='#999999',label='ceiling-bound')],fontsize=8)
fig.savefig(R/'figuras/modelos_q.pdf');fig.savefig(R/'figuras/modelos_q.png',dpi=180);plt.close(fig)
stations=pd.read_csv(ROOT/'domain_stations_q_ampliado.csv').set_index('station')
fig,axs=plt.subplots(1,2,figsize=(9,4),constrained_layout=True)
vlim=max(abs(pd.read_csv(ROOT/'output_V2'/c/'data/sgc_1000_q_station_terms.csv').station_term_s).max() for c in ('sgc_q_mascara_rayos','sgc_q_mascara_rayos_damp400'))
for ax,c,label in zip(axs,('sgc_q_mascara_rayos','sgc_q_mascara_rayos_damp400'),('original damping','adopted damping')):
 t=pd.read_csv(ROOT/'output_V2'/c/'data/sgc_1000_q_station_terms.csv').set_index('station');t=t.join(stations[['lon','lat']]);vals=t.station_term_s-t.station_term_s.mean()
 sc=ax.scatter(t.lon,t.lat,c=vals,cmap='RdBu_r',vmin=-vlim,vmax=vlim,marker='^',s=55,edgecolors='black',linewidths=.3)
 ax.set_title(label);ax.set_xlabel('longitude (°)');ax.set_aspect(1/np.cos(np.radians(7)));ax.grid(alpha=.15)
axs[0].set_ylabel('latitude (°)');fig.colorbar(sc,ax=axs,label='station term relative to network mean (s)',shrink=.8)
fig.savefig(R/'figuras/terminos_estacion.pdf');fig.savefig(R/'figuras/terminos_estacion.png',dpi=180);plt.close(fig)
print(json.dumps(extra,ensure_ascii=True,indent=2),flush=True)
