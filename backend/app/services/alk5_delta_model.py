"""Experimental reference-anchored delta regression; not a deployed activity model."""
from collections import Counter
import numpy as np


def molecular_features(smiles):
    from rdkit import Chem, DataStructs
    from rdkit.Chem import Crippen, Descriptors, QED, rdFingerprintGenerator, rdMolDescriptors
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    features = []
    for text in smiles:
        mol = Chem.MolFromSmiles(text)
        if mol is None:
            raise ValueError('Invalid molecular structure')
        bits = np.zeros(2048, dtype=np.float64)
        DataStructs.ConvertToNumpyArray(generator.GetFingerprint(mol), bits)
        descriptors = [Descriptors.MolWt(mol)/100, Crippen.MolLogP(mol)/2,
            rdMolDescriptors.CalcTPSA(mol)/50, rdMolDescriptors.CalcNumHBA(mol)/3,
            rdMolDescriptors.CalcNumHBD(mol)/2, rdMolDescriptors.CalcNumRotatableBonds(mol)/3,
            QED.qed(mol)/.2]
        features.append(np.r_[bits, descriptors])
    return np.asarray(features)


class ReferenceDeltaModel:
    """Learn within-assay changes, weighting source groups equally.

    Supplied anchor potency must itself be measured; its value never enters fit
    unless explicitly included in training pairs. Support checks are geometric
    heuristics on edit fingerprints, not calibrated uncertainty intervals.
    """
    def fit(self, features, potency, pairs, allowed_indices):
        from sklearn.linear_model import Ridge
        allowed = set(allowed_indices)
        if not pairs:
            raise ValueError('No eligible measured training pairs')
        if any(i not in allowed or j not in allowed for i, j, group in pairs):
            raise ValueError('Training pair crosses allowed training partition')
        differences = np.asarray([features[j]-features[i] for i, j, group in pairs])
        outcomes = np.asarray([potency[j]-potency[i] for i, j, group in pairs])
        self.groups = np.asarray([group for i, j, group in pairs])
        counts = Counter(self.groups.tolist())
        weights = np.asarray([len(pairs)/(len(counts)*counts[g]) for g in self.groups])
        self.model = Ridge(alpha=10, fit_intercept=False, solver='lsqr', tol=1e-5)
        self.model.fit(differences, outcomes, sample_weight=weights)
        bits = differences[:, :2048]
        lengths = np.linalg.norm(bits, axis=1)
        self.unit_edits = bits/np.maximum(lengths[:, None], 1e-12)
        return self

    def predict_change(self, anchor_features, query_features):
        difference = np.asarray(query_features)-np.asarray(anchor_features)
        delta = float(self.model.predict(difference[None])[0])
        bit_change = difference[:2048]
        norm = float(np.linalg.norm(bit_change))
        cosine = np.abs(self.unit_edits@(bit_change/max(norm, 1e-12)))
        support_groups = sorted(set(self.groups[cosine >= .3].tolist())) if norm else []
        a = np.asarray(anchor_features)[:2048]; b = np.asarray(query_features)[:2048]
        union = float(np.maximum(a,b).sum())
        similarity = float(np.minimum(a,b).sum()/union) if union else 1.
        reasons = []
        if similarity < .65:
            reasons.append('anchor_similarity_below_0.65')
        if len(support_groups) < 3:
            reasons.append('edit_supported_by_fewer_than_3_training_source_groups')
        return dict(delta_pic50=delta, anchor_similarity=similarity,
            supporting_training_groups=support_groups,
            max_edit_cosine=float(cosine.max()) if len(cosine) else 0.,
            abstain=bool(reasons), abstention_reasons=reasons,
            uncertainty_status='support_heuristic_not_calibrated_interval')
