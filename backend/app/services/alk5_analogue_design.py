"""Bounded, deterministic peripheral BRICS substitutions, not an activity model."""
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path

ASSET = Path(__file__).resolve().parents[1]/'data/alk5_design_references.json'


def design_analogues(count=30, asset_path=ASSET):
    if not isinstance(count, int) or isinstance(count, bool) or not 1 <= count <= 100:
        raise ValueError('count must be an integer between 1 and 100')
    from rdkit import Chem, DataStructs
    from rdkit.Chem import BRICS, Crippen, Descriptors, QED, rdFingerprintGenerator, rdMolDescriptors
    from rdkit.Chem.MolStandardize import rdMolStandardize
    from rdkit.Chem.Scaffolds import MurckoScaffold

    raw = Path(asset_path).read_bytes()
    asset = json.loads(raw)
    refs = asset['references']
    if not refs or not asset['excluded_normalized_parents']:
        raise ValueError('Missing curated references or identity exclusion archive')
    enum = rdMolStandardize.TautomerEnumerator()
    uncharger = rdMolStandardize.Uncharger()
    def identity(m):
        m = enum.Canonicalize(uncharger.uncharge(rdMolStandardize.FragmentParent(rdMolStandardize.Cleanup(m))))
        Chem.RemoveStereochemistry(m)
        return Chem.MolToSmiles(m)
    def scaffold(m):
        return Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(m))
    def props(m):
        return dict(mw=Descriptors.MolWt(m), logp=Crippen.MolLogP(m), qed=QED.qed(m), tpsa=rdMolDescriptors.CalcTPSA(m), heavy_atoms=m.GetNumHeavyAtoms())
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    mols, fps, scaffold_mols = [], [], []
    for r in refs:
        m = Chem.MolFromSmiles(r['smiles'])
        if m is None or any(str(s.specified) == 'Unspecified' for s in Chem.FindPotentialStereo(m)):
            raise ValueError('Invalid or stereochemically unresolved curated reference')
        core = MurckoScaffold.GetScaffoldForMol(m)
        if core.GetNumHeavyAtoms() < 10:
            raise ValueError('Reference scaffold too small')
        mols.append(m)
        fps.append(generator.GetFingerprint(m))
        scaffold_mols.append(core)
    sides, cores = {}, []
    for i, m in enumerate(mols):
        for bond in sorted(BRICS.FindBRICSBonds(m)):
            broken = BRICS.BreakBRICSBonds(m, bonds=[bond])
            parts = Chem.GetMolFrags(broken, asMols=True)
            if len(parts) != 2:
                continue
            for part in parts:
                dummy_count = sum(a.GetAtomicNum() == 0 for a in part.GetAtoms())
                if dummy_count != 1:
                    continue
                real_atoms = sum(a.GetAtomicNum() > 1 for a in part.GetAtoms())
                text = Chem.MolToSmiles(part)
                if real_atoms <= 12 and not any(a.IsInRing() for a in part.GetAtoms()):
                    sides.setdefault(text, set()).add(i)
                if refs[i].get('design_parent', True) and part.HasSubstructMatch(scaffold_mols[i]):
                    cores.append((i, text, bond[0]))
    # Stable caps bound CPU usage; identical requests have identical output.
    side_items = sorted(sides.items())[:96]
    core_items = sorted(set(cores))[:96]
    seen = set(asset['excluded_normalized_parents'])
    pools = {i: [] for i in range(len(refs))}
    stats = Counter()
    parent_properties = [props(m) for m in mols]
    for i, core_smi, cut in core_items:
        if len(pools[i]) >= 8:
            continue
        core = Chem.MolFromSmiles(core_smi)
        # Rotate deterministic donor order across parents to avoid giving
        # every parent the same alphabetically first replacement fragment.
        offset = i % len(side_items) if side_items else 0
        for side_smi, donor_indices in side_items[offset:]+side_items[:offset]:
            if stats['fragment_pair_attempts'] >= 5000:
                break
            stats['fragment_pair_attempts'] += 1
            side = Chem.MolFromSmiles(side_smi)
            for product in BRICS.BRICSBuild([side], seeds=[core], maxDepth=0, scrambleReagents=False):
                stats['raw_products'] += 1
                try:
                    Chem.SanitizeMol(product)
                except Exception:
                    stats['invalid_products'] += 1
                    continue
                if any(a.GetAtomicNum() == 0 for a in product.GetAtoms()):
                    continue
                stats['sanitized_complete_products'] += 1
                key = identity(product)
                if key in seen:
                    stats['known_or_duplicate'] += 1
                    continue
                if scaffold(product) != scaffold(mols[i]):
                    stats['scaffold_changed'] += 1
                    continue
                if any(str(s.specified) == 'Unspecified' for s in Chem.FindPotentialStereo(product)):
                    stats['unresolved_stereo'] += 1
                    continue
                fp = generator.GetFingerprint(product)
                sim = DataStructs.TanimotoSimilarity(fps[i], fp)
                pp = props(product)
                if not (.5 <= sim <= .9 and pp['mw'] <= 550 and -2 <= pp['logp'] <= 6 and pp['heavy_atoms'] <= mols[i].GetNumHeavyAtoms()+8):
                    stats['design_boundary_reject'] += 1
                    continue
                seen.add(key)
                provenance = dict(parent_id=refs[i]['id'], parent_smiles=refs[i]['smiles'], parent_documents=refs[i]['documents'], parent_activity_ids=refs[i]['activity_ids'], parent_scaffold=scaffold(mols[i]), retained_fragment=core_smi, replacement_fragment=side_smi, donor_ids=[refs[j]['id'] for j in sorted(donor_indices)], parent_cut_atom_indices=list(cut), transformation='single peripheral BRICS fragment substitution')
                pools[i].append(dict(smiles=Chem.MolToSmiles(product), normalized_parent=key, parent_similarity=sim, max_reference_similarity=max(DataStructs.BulkTanimotoSimilarity(fp, fps)), properties=pp, property_changes={k: pp[k]-parent_properties[i][k] for k in pp}, provenance=provenance, activity_status='unmeasured', safety_status='not_evaluated', synthesis_status='route_not_established'))
                if len(pools[i]) >= 8:
                    break
            if len(pools[i]) >= 8:
                break
    # Round-robin across parents; output position is not activity ranking.
    candidates = []
    for offset in range(8):
        for i in pools:
            if offset < len(pools[i]) and len(candidates) < count:
                candidates.append(dict(candidate_id=f'ALK5_design_{len(candidates)+1:03d}', **pools[i][offset]))
    candidate_fps = [generator.GetFingerprint(Chem.MolFromSmiles(r['smiles'])) for r in candidates]
    pair_sims = [DataStructs.TanimotoSimilarity(candidate_fps[i], candidate_fps[j]) for i in range(len(candidate_fps)) for j in range(i)]
    return dict(target='ALK5', design_method='reference_anchored_brics_v1', requested=count, returned=len(candidates), candidates=candidates,
        metrics=dict(reference_count=len(refs), design_parent_count=sum(r.get('design_parent', True) for r in refs), parent_families_represented=len({r['provenance']['parent_id'] for r in candidates}), unique_scaffolds=len({r['provenance']['parent_scaffold'] for r in candidates}), replacement_fragments_used=len({r['provenance']['replacement_fragment'] for r in candidates}), enumerated_unique_designs=sum(len(p) for p in pools.values()), mean_pairwise_tanimoto=sum(pair_sims)/len(pair_sims) if pair_sims else None, counts=dict(stats)),
        meta=dict(reference_asset_sha256=sha256(raw).hexdigest(), reference_provenance=asset['provenance'], status='design_hypotheses_only', activity_validated=False, safety_validated=False, order='round_robin_by_parent_not_activity', exclusion_scope='archived normalized structures only; not global or patent novelty', admet='not_run', docking='not_run', model_training='none; rule-based fragment enumeration', limitations=['Same scaffold and similarity do not establish activity or selectivity', 'BRICS compatibility does not establish an executable synthetic route', 'Finite library; requested count is not guaranteed']))
