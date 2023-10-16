# Shape Detection and Mediation Analysis using Semi-parametric Shape-Restricted Regression Spline
## Qing Yin

This R package is designed according to my doctoral dissertation "Shape Detection and Mediation Analysis using Semi-parametric Shape-Restricted Regression Spline with Applications".

In the first part of the dissertation, a shape detection method based on regression splines is developed. The proposed method can help researchers select the most suitable shape to describe their data among increasing, decreasing, convex and concave shapes. Specifically, we develop a technique based on mixed effects regression spline to analyze hormonal data, but the method is general enough to be applied to other similar problems.

In the second part of the dissertation, we develop a method to analytically estimate the direct and indirect effects when we have some prior knowledge on the relationship between the mediator and the outcome (increasing, decreasing, convex or concave). In order to make suitable inferences, the asymptotic confidence intervals of those effects are obtained via delta method.

To install the package, please follow the following steps:
```{r, eval = F}
library(devtools)
install_github("qy2166/SRMERS")
```

To load the package into your local enviroment, please follow the following step:
```{r, eval = F}
library(SRMERS)
```

Before conducting the analysis, please follow the following suggestions regarding your data set:

    1. The data set should be complete.
    2. You can use the function model.matrix() to transform factors into dummy variables
    3. You can use the function data.frame() to make the data set into a data frame

To perform the shape detection using a fixed effect model, you may use the function FERS(). Please see the example:
```{r, eval = F}
shape <- FERS(y = "ySim",
              xMain = "hormone",
              xConf = c("age", "invwt", "race2", "race3", "race4", "race5",
                        "season2", "season3", "season4", "smoking1", "ovum1", "diabetes1"),
              dataset = data.sim.fixed,
              nBasis = 5,
              nIter = 1000)
```

To perform the shape detection using a mixed effects model, you may use the function MERS(). Please see the example:
```{r, eval = F}
shape <- MERS(y = "ySim",
              xMain = "hormone",
              xConf = c("age", "invwt", "race2", "race3", "race4", "race5",
                        "season2", "season3", "season4", "smoking1", "ovum1", "diabetes1"),
              xRand = "cluster",
              dataset = data.sim.mixed,
              nBasis = 5,
              nIter = 1000)
```

To perform the mediation analysis using a fixed effect model, you may use the function SRSplineMed(). Please see the example:
```{r, eval = F}
medModel <- SRSplineMed(data = data.sim.med, nBasis = 5, 
                        exposure = "pesticide1",
                        mediator = "hormone",
                        outcome = "ySim",
                        confounderVec = c("age", "invwt", "race2", "race3", "race4", "race5",
                                          "season2", "season3", "season4", "smoking1", "ovum1",
                                          "diabetes1"),
                        shapeExp = "concave", shapeNonExp = "increasing",
                        mValue = 0.15, varAsymp = T)
```

You may also use the functions LRMed (binary exposure) or LRmed2 (continuous exposure) to perform the mediation analysis using the VanderWeele's approach. Please see the example:
```{r, eval = F}
medModel <- LRMed(data = data.sim.med, 
                  exposure = "pesticide1", 
                  mediator = "hormone", 
                  outcome = "ySim", 
                  confounderVec = c("age", "invwt", "race2", "race3", "race4", "race5",
                                    "season2", "season3", "season4", "smoking1", "ovum1",
                                    "diabetes1"), 
                  mValue = 0.15)
```

Acknowledgement:

Thank you for the support and help from my advisors, Dr. Shyamal D. Peddada, Dr. Jennifer J. Adibi, and Dr. Jong H. Jeong, at University of Pittsburgh in developing this R package.

References:

    1. Yin Q, Xun X, Peddada SD, Adibi JJ. Shape Detection Using Semi-Parametric Shape-Restricted Mixed Effects Regression Spline with Applications. Sankhya B 2021; 83(S1): 65-85.
    2. VanderWeele TJ. Explanation In Causal Inference : Methods for Mediation and Interaction. Oxford University Press, New York. 2015.
    3. Yin Q, Jeong JH, Qin X, Peddada SD, Adibi JJ. Mediation Analysis using Semi-parametric Shape-Restricted Regression with Applications. arXiv:2310.09185 [stat.ME] 2023.
