"""Render numerical audit outputs into TeX; never modify archived campaigns."""
import json
from pathlib import Path
R=Path(__file__).resolve().parent
j=json.loads((R/'segunda_revision/diagnosticos.json').read_text())
s=j['bootstrap']['campaigns']; out=[]
def macro(name,value):out.append('\\newcommand{\\'+name+'}{'+str(value)+'}')
for key,v in s.items():
 for field,name in [('E_RMS','E'),('structure_ratio','Structure'),('station_ratio','Station'),('component_corr','ComponentCorr'),('prediction_data_corr','PredictionCorr')]:
  macro('audit'+key+name,f'{v[field]:.4f}')
 ci=v['conditional_event_bootstrap_95CI'];macro('audit'+key+'CI',f'{ci[0]:.2f}--{ci[1]:.2f}')
 for field,name in [('initial_residual_rms_s','InitialRms'),('data_rms_s','DataRms')]:macro('audit'+key+name,f'{v[field]:.5f}')
macro('auditBootstrapReplicates',j['bootstrap']['replicates']);macro('auditBootstrapSeed',j['bootstrap']['seed'])
for key,name in [('RayDomain_minus_Expanded','DomainGain'),('Adopted_minus_RayDomain','DampingChange')]:
 v=j['bootstrap']['paired_contrasts'][key];macro('audit'+name,f"{v['difference_percentage_points']:.2f}")
 ci=v['paired_event_bootstrap_95CI'];macro('audit'+name+'CI',f'{ci[0]:.2f} to {ci[1]:.2f}')
for v in j['controlled_experiment']['fits']:
 macro('audit'+v['velocity']+v['mask']+'E',f"{v['E_RMS']:.4f}")
 macro('audit'+v['velocity']+v['mask']+'Rms',f"{v['residual_rms_s']:.5f}")
fits={(v['velocity'],v['mask']):v for v in j['controlled_experiment']['fits']}
macro('auditMaskGainFixedVp',f"{fits['InheritedVp','RayMask']['E_RMS']-fits['InheritedVp','InheritedMask']['E_RMS']:.4f}")
macro('auditVelocityGainRayMask',f"{fits['UpdatedVp','RayMask']['E_RMS']-fits['InheritedVp','RayMask']['E_RMS']:.4f}")
d=json.loads((R/'revision_cientifica/diagnosticos.json').read_text())
def find(obj,key):
 if isinstance(obj,dict):
  if key in obj:return obj[key]
  for v in obj.values():
   found=find(v,key)
   if found is not None:return found
 return None
for key,name in [('clipped_q_rms_300_s','ClippedLowE'),('clipped_q_rms_1500_s','ClippedHighE')]:
 macro('audit'+name,f"{100*(1-find(d,key)/s['Adopted']['data_rms_s']):.2f}")
macro('auditLowerBoundNodes',j['bounds']['Adopted']['at_lower_bound'])
macro('auditLowerBoundPercent',f"{j['bounds']['Adopted']['lower_bound_fraction_percent']:.1f}")
(R/'cifras_segunda_revision.tex').write_text('% Generated from segunda_revision/diagnosticos.json\n'+'\n'.join(out)+'\n',encoding='utf-8')
print('Generated cifras_segunda_revision.tex')
