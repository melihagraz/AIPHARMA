"""Experimental assay-offset Gaussian process conditioned on a measured reference."""
import numpy as np
from scipy.linalg import cho_factor, cho_solve


def tanimoto(first, second):
    first=np.asarray(first,dtype=float); second=np.asarray(second,dtype=float)
    product=first@second.T
    denominator=first.sum(axis=1)[:,None]+second.sum(axis=1)[None,:]-product
    return product/np.maximum(denominator,1e-12)


class AssayReferenceGP:
    """Fixed kernel and noise; intervals are model-based until empirically checked.

    The assay offset represents nuisance assay variability, not a biological
    mechanism. The new assay must not appear in fit; the supplied reference
    measurement enters posterior conditioning only.
    """
    noise_variance=.04
    assay_variance=1.

    def fit(self,bits,potency,assays):
        self.bits=np.asarray(bits,dtype=float)
        self.assays=np.asarray(assays)
        self.mean=float(np.median(potency))
        kernel=tanimoto(self.bits,self.bits)+self.assay_variance*(self.assays[:,None]==self.assays[None,:])
        kernel.flat[::len(kernel)+1]+=self.noise_variance
        self.factor=cho_factor(kernel,lower=True,check_finite=True)
        self.alpha=cho_solve(self.factor,np.asarray(potency)-self.mean)
        return self

    def predict_with_reference(self,query,reference,reference_pic50,new_assay):
        if new_assay in self.assays:
            raise ValueError('This evaluation requires a held-out assay')
        query=np.asarray(query,dtype=float)[None,:]
        reference=np.asarray(reference,dtype=float)[None,:]
        kq=tanimoto(query,self.bits)[0]; ka=tanimoto(reference,self.bits)[0]
        sq=cho_solve(self.factor,kq); sa=cho_solve(self.factor,ka)
        mean_q=self.mean+float(kq@self.alpha); mean_a=self.mean+float(ka@self.alpha)
        similarity=float(tanimoto(query,reference)[0,0])
        cross=similarity+self.assay_variance-float(kq@sa)
        var_a=1+self.assay_variance-float(ka@sa)+self.noise_variance
        var_q=1+self.assay_variance-float(kq@sq)
        prediction=mean_q+cross/var_a*(reference_pic50-mean_a)
        variance=var_q-cross**2/var_a+self.noise_variance
        if variance < -1e-8:
            raise ValueError('Invalid negative conditional variance')
        sd=float(np.sqrt(max(variance,0)))
        reasons=[]
        if similarity<.65: reasons.append('reference_similarity_below_0.65')
        if sd>.7: reasons.append('model_standard_deviation_above_0.7_pic50')
        return dict(prediction=float(prediction),standard_deviation=sd,
            interval_95=[float(prediction-1.96*sd),float(prediction+1.96*sd)],
            reference_similarity=similarity,abstain=bool(reasons),abstention_reasons=reasons,
            interval_status='model_based_not_independently_calibrated',activity_validated=False)
