import torch
import torch.nn as nn
import zuko

class SingleFlowMorpher(nn.Module):
    def __init__(self, num_features=4, num_condition=4, num_transforms=4, bins=10, hidden_features=[64, 64]):
        """
        num_features: Number of variables to transform (number of informative features)
        num_condition:  Number of conditioning variables ( IsData + ancillary features)
        hidden_features: structure of the pperceptron (in this case 2-layered multiperceptron with each layer having 64 units)
        """
        super().__init__()
        
        # zuko NSF uses Masked Autoregressive Flows (MADE) with monotonic rational quadratic splines (ref: https://zuko.readthedocs.io/0.3.2/api/zuko.flows.spline.html)
        self.flow = zuko.flows.NSF(
            features=num_features,
            context=num_condition,
            transforms=num_transforms,
            bins=bins,
            hidden_features=hidden_features
        )

    def forward(self, x, context):
        """
        Computing NLL (Loss function) during training
        """
        # zuko's flow(condition) returns a distribution object
        # .log_prob(x) calculates -(log pz(f(x)) + log|det J|)
        return -self.flow(context).log_prob(x).mean()

    def transform_to_base(self, x, context):
        """
        Forward Pass (f): Maps input feature x to Gaussian base space z.
        """
        # forward bijection f(x)
        return self.flow(context).transform(x)

    def transform_from_base(self, z, context):
        """
        Inverse Pass (f^-1): Maps Gaussian (base) point z back to feature space x.
        """
        # inverse bijection f^{-1}(z)
        return self.flow(context).transform.inv(z)
