# R Notebook


```r
library(dplyr)
```

```
## Warning: package 'dplyr' was built under R version 3.4.2
```

```
## 
## Attaching package: 'dplyr'
```

```
## The following objects are masked from 'package:stats':
## 
##     filter, lag
```

```
## The following objects are masked from 'package:base':
## 
##     intersect, setdiff, setequal, union
```

```r
library(ggplot2)

MB231 <- read.csv("./data/MB231.csv", stringsAsFactors = FALSE) %>%
  tidyr::separate(X, into = c("matrix", "drug"), sep = " ") %>%
  mutate(matrix = as.factor(matrix)) %>%
  dplyr::select(prolifTwo, CREB, EGFR, NFkB, p38, AKT, JNK, MEK1, ERK1.2, matrix, drug)

rescaleLMsum <- function(dataIn) {
  dataIn <- select(dataIn, -matrix, -drug) %>% scale(center = TRUE, scale = TRUE) %>% data.frame

  # Applying the linear model
  model <- lm(prolifTwo ~ . + 0, data = dataIn)
  
  # Summarize the model
  print(summary(model))
}

remove3Dmatrix <- function(dataIn) {
  return(dataIn %>% filter(matrix != "sph") %>% filter(matrix != "single") %>% filter(matrix != "serum"))
}

removeOtherDrug <- function(dataIn) {
  return(dataIn %>% filter(drug != "lap") %>% filter(drug != "tems"))
}
```

The full model:


```r
# Scaling the phospho-measurements by z-score
MB231 %>% rescaleLMsum()
```

```
## 
## Call:
## lm(formula = prolifTwo ~ . + 0, data = dataIn)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -1.5457 -0.4924 -0.1094  0.4007  3.1405 
## 
## Coefficients:
##        Estimate Std. Error t value Pr(>|t|)    
## CREB    0.08839    0.19838   0.446   0.6573    
## EGFR   -0.52229    0.22045  -2.369   0.0205 *  
## NFkB    0.10901    0.18457   0.591   0.5566    
## p38     0.08970    0.21109   0.425   0.6722    
## AKT    -0.38011    0.15950  -2.383   0.0198 *  
## JNK     0.28277    0.16948   1.669   0.0996 .  
## MEK1    0.94257    0.21496   4.385 3.89e-05 ***
## ERK1.2 -0.31953    0.18466  -1.730   0.0878 .  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.8413 on 72 degrees of freedom
## Multiple R-squared:  0.3549,	Adjusted R-squared:  0.2832 
## F-statistic: 4.952 on 8 and 72 DF,  p-value: 6.646e-05
```

The model without lap or tems:


```r
# Scaling the phospho-measurements by z-score
MB231 %>% removeOtherDrug() %>% rescaleLMsum()
```

```
## 
## Call:
## lm(formula = prolifTwo ~ . + 0, data = dataIn)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -1.8109 -0.6387  0.3275  0.6592  1.2066 
## 
## Coefficients:
##         Estimate Std. Error t value Pr(>|t|)  
## CREB    0.447504   0.407460   1.098   0.2803  
## EGFR   -0.911214   0.444173  -2.051   0.0485 *
## NFkB    0.210109   0.354375   0.593   0.5574  
## p38     0.175151   0.250221   0.700   0.4890  
## AKT    -0.554166   0.290082  -1.910   0.0651 .
## JNK     0.210388   0.214338   0.982   0.3337  
## MEK1    0.697828   0.338423   2.062   0.0474 *
## ERK1.2 -0.001727   0.400379  -0.004   0.9966  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.9157 on 32 degrees of freedom
## Multiple R-squared:  0.312,	Adjusted R-squared:   0.14 
## F-statistic: 1.814 on 8 and 32 DF,  p-value: 0.1109
```

The model without sph or single matrices:


```r
# Scaling the phospho-measurements by z-score
MB231 %>% remove3Dmatrix() %>% rescaleLMsum()
```

```
## 
## Call:
## lm(formula = prolifTwo ~ . + 0, data = dataIn)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -1.3178 -0.3737 -0.2106  0.2604  1.8775 
## 
## Coefficients:
##        Estimate Std. Error t value Pr(>|t|)   
## CREB    0.29255    0.31932   0.916  0.36869   
## EGFR   -0.53961    0.32455  -1.663  0.10939   
## NFkB   -0.48681    0.33338  -1.460  0.15720   
## p38     0.14134    0.37038   0.382  0.70612   
## AKT     0.11846    0.28198   0.420  0.67815   
## JNK     0.07434    0.26419   0.281  0.78082   
## MEK1    0.91314    0.30805   2.964  0.00676 **
## ERK1.2 -0.53583    0.26986  -1.986  0.05862 . 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.8476 on 24 degrees of freedom
## Multiple R-squared:  0.4438,	Adjusted R-squared:  0.2584 
## F-statistic: 2.394 on 8 and 24 DF,  p-value: 0.04691
```

The model without single sph matrices, lab or tems:


```r
# Scaling the phospho-measurements by z-score
MB231 %>% removeOtherDrug() %>% remove3Dmatrix() %>% rescaleLMsum()
```

```
## 
## Call:
## lm(formula = prolifTwo ~ . + 0, data = dataIn)
## 
## Residuals:
##      Min       1Q   Median       3Q      Max 
## -0.55793 -0.11577 -0.02999  0.17194  0.36331 
## 
## Coefficients:
##        Estimate Std. Error t value Pr(>|t|)   
## CREB    0.80294    0.26714   3.006  0.01693 * 
## EGFR   -1.00300    0.23987  -4.181  0.00307 **
## NFkB   -0.61279    0.33590  -1.824  0.10556   
## p38     0.28211    0.18158   1.554  0.15886   
## AKT     0.31647    0.19450   1.627  0.14237   
## JNK     0.08242    0.10588   0.778  0.45870   
## MEK1    0.15466    0.30698   0.504  0.62799   
## ERK1.2  0.29008    0.35970   0.806  0.44330   
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.3228 on 8 degrees of freedom
## Multiple R-squared:  0.9444,	Adjusted R-squared:  0.8888 
## F-statistic: 16.99 on 8 and 8 DF,  p-value: 0.0002916
```



