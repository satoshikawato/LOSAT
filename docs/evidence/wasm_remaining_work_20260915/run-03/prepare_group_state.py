# NCBI reference: c++/src/algo/blast/core/link_hsps.c:603-688,827-895,955-980
# H2_index=next_larger; H_hsp_link=H2; H->hsp->evalue=evalue;
# Instrument separate frozen sources, never production or timing artifacts.
from pathlib import Path
import shutil, json, hashlib, re, argparse
e=Path(__file__).resolve().parent
p=argparse.ArgumentParser()
p.add_argument('--versions',nargs='+',default=['baseline','I1'])
a=p.parse_args()
template=(e/'group_state_template.rs').read_text()
for version in a.versions:
    dest=e/'work'/(version+'state-source')
    source=e/'work'/(version+'-source')
    if dest.exists():
        for p in (source/'LOSAT/src').rglob('*.rs'):
            assert p.read_bytes()==(dest/p.relative_to(source)).read_bytes(),p
    else:
        shutil.copytree(source,dest)
    p=dest/'LOSAT/src/algorithm/tblastx/sum_stats_linking/linking.rs'
    s=p.read_text()
    assert s.count('#[cfg(test)]\nmod tests {')==1
    s=s.replace('#[cfg(test)]\nmod tests {','mod tests {')
    s=s.replace('    fn mock_hit(','    pub(super) fn mock_hit(')
    s+='\n'+template
    needle='    while remaining > 0 {'
    assert s.count(needle)==1
    s=s.replace(needle,needle+'\n        println!("STATE ENTER {:?}", (selection_round,remaining,first_pass,path_changed,state_index(active_head)));\n        state_snapshot("ENTER",pool_lh_helpers,pool_hsp_links,&group_hits,&maximum_tree);')
    needle='        selection_round += 1;'
    assert s.count(needle)==1
    s=s.replace(needle,'        println!("STATE EXIT {:?}", (selection_round,remaining,first_pass,path_changed,state_index(active_head),best,best_sum,use_current_max,prob.map(f64::to_bits),ordering,evalue.to_bits(),linked_set));\n        state_snapshot("EXIT",pool_lh_helpers,pool_hsp_links,&group_hits,&maximum_tree);\n'+needle)
    # NCBI link_hsps.c:827-863: record the same H/H2 boundary when the
    # unchanged visit body is a module macro. Pass the owner explicitly so
    # macro hygiene cannot bind diagnostic state to a different caller.
    if 'macro_rules! visit_large_gap_predecessor' in s:
        needle='$next_cursor:expr;\n     $h_sum:ident'
        assert s.count(needle)==1
        s=s.replace(needle,'$next_cursor:expr;\n     $state_owner:expr; $h_sum:ident')
        pattern=r'(visit_large_gap_predecessor!\((?:prefix|indexed),[^;]+;\s*)(h_sum,)'
        assert len(re.findall(pattern,s))==3
        s=re.sub(pattern,r'\1pool_lh_helpers[h_lh_idx].hsp_idx; \2',s)
        needle='let b0 = sum <= $h_sum;'
        assert s.count(needle)==1
        s=s.replace(needle,needle+'\nprintln!("STATE VISIT {} {} {} {} {} {} {:?} {}",$state_owner,$current_idx,next_larger,b0,$h_sum,$h_num,state_index($h_link),$h_xsum.to_bits());')
    else:
        needle='let b0 = sum <= h_sum;'
        assert s.count(needle)==1
        owner='i' if version=='baseline' else 'pool_lh_helpers[h_lh_idx].hsp_idx'
        s=s.replace(needle,needle+'\nprintln!("STATE VISIT {} {} {} {} {} {} {:?} {}",'+owner+',current_idx,next_larger,b0,h_sum,h_num,state_index(h_link),h_xsum.to_bits());')
    needle='let new_sum = h_sum + (i_score - cutoff_big);'
    assert s.count(needle)==1
    s=s.replace(needle,'println!("STATE CHOICE {} {} {} {:?} {}",i,h_sum,h_num,state_index(h_link),h_xsum.to_bits());\n'+needle)
    # Deliberately force the existing per-HSP trace branch in these diagnostic
    # copies only; off/on exercise identical algorithm state at every boundary.
    pattern=r'let is_target_hsp = if debug_chaining \{\s*target_hsp_idx == Some\(i\)\s*\} else \{\s*false\s*\};'
    assert len(re.findall(pattern,s))==2
    s=re.sub(pattern,'let is_target_hsp = std::env::var_os("LOSAT_STATE_TRACE").is_some();',s)
    p.write_text(s)
    p=dest/'LOSAT/src/algorithm/tblastx/sum_stats_linking/mod.rs'
    p.write_text(p.read_text()+'\n// NCBI link_hsps.c:603-982: while(number_of_hsps>0) { ... }\n// Diagnostic-only entry for frozen group transcripts.\npub use linking::remaining_state_transcript;\n')
    (dest/'LOSAT/src/main.rs').write_text('// NCBI link_hsps.c:603-982: while(number_of_hsps>0) { ... }\n// Diagnostic-only entry, not a LOSAT runtime feature.\nfn main() { LOSAT::algorithm::tblastx::sum_stats_linking::remaining_state_transcript(); }\n')
    manifest={str(p.relative_to(dest)):hashlib.sha256(p.read_bytes()).hexdigest() for p in dest.rglob('*') if p.is_file()}
    (e/(version+'state-source-manifest.json')).write_text(json.dumps(manifest,indent=2)+'\n')
    print(version,'state source prepared',flush=True)
