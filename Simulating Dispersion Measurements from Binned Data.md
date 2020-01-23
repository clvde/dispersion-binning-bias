### Seed Dispersion Data Collection Simulation

this script goes through data simulation for seed dispersal and measures the bias on dispersion measurements where the angle is binned.

Courtney Van Den Elzen, January 2020

#### Custom functions


```r
ptCentDist <- function(xpoint,ypoint,xcent,ycent) {
sqrt(abs(xpoint-xcent)^2 + abs(ypoint-ycent)^2)
}
```

#### Define some constants


```r
n_seeds <- 20 # number of seed positions to simulate
shp <- 3
rt <- 1
bin_size <- 22.5 # bin size for the slices of the angle pie. 22.5 degree bins
samps <- 1000 # number of samples to do the calculation on
```

Define the output vecots for the loop


```r
mean_disp_cts_vec <- rep(NA, samps)
mean_disp_binned_vec  <- rep(NA, samps)
```

#### Turn this into a loop and output distributions of each


```r
for (s in 1:samps){
  #' distances drawn from a normal 
  dists <- rgamma(n_seeds, shp, rt)
  #' Angles drawn from a uniform
  angles <- runif(n_seeds, min = 0, max = 360)
  #' create angle bins and bin angle values (middle of the slice, so bin 1 = mean(0,22.5) deg = 11.25 deg)
  angle_bins <- floor(1 + angles / 22.5)
  bin_mid_angle <- angle_bins * 22.5 - 22.5/2

  #' Create a data frame of distances and angles      
  df <- data.frame(dists, angles, angle_bins, bin_mid_angle)

  #' Convert to cartesian coordinates
  #' Distance values in rad
  r <- df$dists

  #' Continuous version - not binning angles
  df$dist_x_cts <- r*cos(df$angles)
  df$dist_y_cts <- r*sin(df$angles)
  centroid_cts <- c(mean(df$dist_x_cts), mean(df$dist_y_cts))

  #' Binned version - using midpoint values of angle bins
  df$dist_x_binned <- r*cos(df$bin_mid_angle)
  df$dist_y_binned <- r*sin(df$bin_mid_angle)
  centroid_binned <- c(mean(df$dist_x_binned), mean(df$dist_y_binned))

  #' Calculate the dispersion of points
  disp_vec_cts <- ptCentDist(df$dist_x_cts, centroid_cts[1], df$dist_y_cts, centroid_cts[2])
  disp_vec_binned <- ptCentDist(df$dist_x_binned, centroid_binned[1], df$dist_y_binned, centroid_binned[2])
  
  #' Put mean dispersion values into vectors
  mean_disp_cts_vec[s] <- mean(disp_vec_cts)
  mean_disp_binned_vec[s] <- mean(disp_vec_binned)
}
```

#### Results

Linear model - what is the relationship between binned and cts average dispersion values?


```r
summary(lm(mean_disp_binned_vec~mean_disp_cts_vec))
```

```
## 
## Call:
## lm(formula = mean_disp_binned_vec ~ mean_disp_cts_vec)
## 
## Residuals:
##      Min       1Q   Median       3Q      Max 
## -1.42599 -0.36010 -0.04801  0.30456  2.06801 
## 
## Coefficients:
##                   Estimate Std. Error t value Pr(>|t|)    
## (Intercept)        1.54532    0.09422    16.4   <2e-16 ***
## mean_disp_cts_vec  0.54600    0.03212    17.0   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.5207 on 998 degrees of freedom
## Multiple R-squared:  0.2245,	Adjusted R-squared:  0.2237 
## F-statistic: 288.9 on 1 and 998 DF,  p-value: < 2.2e-16
```

Correlation between binned and cts average dispersion vectors


```r
cor(mean_disp_binned_vec,mean_disp_cts_vec)
```

```
## [1] 0.4737952
```

Plot the relationship between average binned and cts dispersion values


```r
disp_df <- data.frame(mean_disp_binned_vec, mean_disp_cts_vec)

ggplot(data = disp_df, mapping = aes(x = mean_disp_cts_vec, y = mean_disp_binned_vec)) +
  geom_point() + 
  geom_smooth(method='lm')
```

![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-7-1.png)

