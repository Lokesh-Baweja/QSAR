# Quantitative Structure Activity Relationship (QSAR) in R to understand the molecular properties responsible for mRNA Binding

The data and code for this study was taken from a recent study "Quantitative Structure–Activity Relationship (QSAR) Study Predicts Small-Molecule Binding to RNA Structure" https://pubs.acs.org/doi/10.1021/acs.jmedchem.2c00254. I have compiled a jupyter notebook to reproduce the model building. In the dataset, we have 48 molecules with 400+ decriptors. lnKD (dissociation constant is the dependent variable) 
So, this is a situation where number of data points (N) << (Number of Features or descriptors). Classical example of curse of dimensionality. 
Multiple Linear Regression seems to be a best approach for building model building. What important here is how we can reduce the dimension of the dataset. One of the approaches is to find best features for building predictive models in LASSO (Least Absolute Shrinkage and Selection operator) is an efficient method aka L1 regularization. A more detailed description can be found at https://www.publichealth.columbia.edu/research/population-health-methods/least-absolute-shrinkage-and-selection-operator-lasso.

Library Requirements for performing this analysis

"prospectr"

"glmnet"

## Steps for performing this analysis

1. Load the preprocessed dataset and perform the LASSO analysis to find the  number of features based on the minimum mean-sqaured error plotted against the lambda values. 

2. One of the challenges with this approach is you have to generate different combination of features from the lasso selected features to get the best model exlaining the structure activity relationship. This is not a trivial task. Domain knowledge plays a key role in selecting the appropriate model

2. Using the best model to predict the lnKD values



```R
data_qsar <- read.csv("KD_refine.csv")
```


```R
head(data_qsar)
```


<table class="dataframe">
<caption>A data.frame: 6 × 194</caption>
<thead>
	<tr><th></th><th scope=col>lnKD</th><th scope=col>rgyr</th><th scope=col>glob</th><th scope=col>AM1_dipole</th><th scope=col>AM1_E</th><th scope=col>AM1_HF</th><th scope=col>AM1_HOMO</th><th scope=col>apol</th><th scope=col>ASA.</th><th scope=col>ASA..1</th><th scope=col>⋯</th><th scope=col>vsurf_R</th><th scope=col>vsurf_W2</th><th scope=col>vsurf_W8</th><th scope=col>vsurf_Wp1</th><th scope=col>vsurf_Wp2</th><th scope=col>vsurf_Wp4</th><th scope=col>vsurf_Wp5</th><th scope=col>vsurf_Wp6</th><th scope=col>vsurf_Wp7</th><th scope=col>weinerPath</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>⋯</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;int&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>1</th><td>-16.26969</td><td>4.820693</td><td>0.1465251</td><td>21.54953</td><td>-207757.5</td><td>236.24567</td><td>-17.88922</td><td>90.86413</td><td>405.9668</td><td>122.86533</td><td>⋯</td><td>1.516588</td><td>1579.6133</td><td>13.565070</td><td>1592.945</td><td>1038.1184</td><td>215.3457</td><td>79.81318</td><td>18.71441</td><td> 3.248806</td><td>6111</td></tr>
	<tr><th scope=row>2</th><td>-12.51714</td><td>4.812603</td><td>0.1497168</td><td>19.21120</td><td>-210033.0</td><td>-42.75333</td><td>-16.38121</td><td>89.35196</td><td>388.9543</td><td>134.00852</td><td>⋯</td><td>1.523779</td><td>1571.8186</td><td>14.939443</td><td>1519.623</td><td> 963.3679</td><td>191.6073</td><td>69.07973</td><td>14.05777</td><td> 2.066936</td><td>6111</td></tr>
	<tr><th scope=row>3</th><td>-12.45121</td><td>4.524976</td><td>0.1029352</td><td>25.71803</td><td>-143880.0</td><td>542.27912</td><td>-18.58516</td><td>71.85307</td><td>317.2251</td><td> 65.44562</td><td>⋯</td><td>1.454466</td><td>1049.0908</td><td> 5.502652</td><td>1313.306</td><td> 814.9738</td><td>164.1943</td><td>62.12614</td><td>16.36669</td><td> 4.427319</td><td>2700</td></tr>
	<tr><th scope=row>4</th><td>-11.41974</td><td>4.199010</td><td>0.4194048</td><td>25.76304</td><td>-196863.4</td><td>-60.47815</td><td>-15.34944</td><td>81.90693</td><td>343.7094</td><td>171.51191</td><td>⋯</td><td>1.528474</td><td>1406.1336</td><td>15.357475</td><td>1350.648</td><td> 770.1565</td><td>168.3968</td><td>83.09161</td><td>41.93077</td><td>20.083616</td><td>4989</td></tr>
	<tr><th scope=row>5</th><td>-12.54856</td><td>4.508194</td><td>0.1117372</td><td>24.42179</td><td>-155728.6</td><td>317.78826</td><td>-17.55147</td><td>71.45014</td><td>342.6224</td><td> 87.98411</td><td>⋯</td><td>1.450390</td><td>1313.0281</td><td> 9.607074</td><td>1377.419</td><td> 897.8370</td><td>184.2273</td><td>65.71518</td><td>13.31651</td><td> 1.539569</td><td>2932</td></tr>
	<tr><th scope=row>6</th><td>-12.48165</td><td>4.124427</td><td>0.1990713</td><td>14.22592</td><td>-151720.2</td><td>443.54318</td><td>-18.51158</td><td>79.24447</td><td>216.7050</td><td> 35.61854</td><td>⋯</td><td>1.547151</td><td> 874.4754</td><td> 9.556666</td><td>1196.090</td><td> 671.2119</td><td>135.7255</td><td>52.62777</td><td>13.87811</td><td> 2.658363</td><td>3234</td></tr>
</tbody>
</table>




```R
#Let's create a function for evaluating results 
#The function takes true and predicted as 

evaluate_results <- function(true, predicted, df) {
    SSE <- ((true-predicted)^2)
    SST <- ((true-mean(true))^2)
    R2 <- 1 - SSE/SST
    RMSE <- sqrt(SSE/nrow(df))
            
    data.frame(
        RMSE = RMSE, 
        R_squared = R2)
    
}

```


```R
#remotes::install_github("l-ramirez-lopez/prospectr")
```


```R
library(prospectr)
```


```R
xspace <- data_qsar[,-1]
ks <- kenStone(as.matrix(xspace), k=12, metric = "mahal",pc=0.99, .center =
TRUE, .scale = FALSE)
ks$test
trainid <- ks$test
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>1</li><li>3</li><li>4</li><li>5</li><li>7</li><li>8</li><li>9</li><li>11</li><li>12</li><li>14</li><li>15</li><li>16</li><li>17</li><li>18</li><li>20</li><li>21</li><li>22</li><li>23</li><li>25</li><li>26</li><li>28</li><li>29</li><li>30</li><li>31</li><li>33</li><li>35</li><li>37</li><li>38</li><li>40</li><li>41</li><li>42</li><li>43</li><li>44</li><li>45</li><li>46</li><li>47</li></ol>




```R
#After performing let's divide the data into training and testset

training_set <- data_qsar[trainid,]
test_set <- data_qsar[-trainid, ]


x_train <- as.matrix(training_set[-1])
y_train <- data.matrix(training_set[1])
x_test <- as.matrix(test_set[-1])
y_test <- as.matrix(test_set[1])

```


```R

```


```R
#install.packages("glmnet")
```


```R
library("glmnet")
```


```R
set.seed(1)
lambdas <- 10^seq(2, -6, length = 100)
print (length(lambdas))
# use sv.glmnet to find the best lambda for lasso from 5-fold cv
lasso_reg <- cv.glmnet(x_train, y_train, alpha = 1, lambda = lambdas,
standardize = TRUE, nfolds = 5)
plot(lasso_reg)
# plot the shrinkage graph with multiple lambda values
lasso_model <- glmnet(x_train, y_train, alpha = 1, nlambda =100,standardize =
TRUE)

p1 <- plot(lasso_model,xvar="lambda",label = T, lwd=4,cex.lab=
2,cex.axis=2,xlim = c(-4.5,0.5), ylim=c(-20,20))

p1.lty=2
box(lwd=4)
# chose the lambda with lowest mean-squrared error from cv
lambda_best_lasso <- lasso_reg$lambda.min

# build the lasso regression model using selected descriptors
lasso_model <- glmnet(x_train, y_train, alpha = 1, lambda
=lambda_best_lasso,standardize = TRUE)

# find the non-zero coefficients and their names
lasso.coef <- predict(lasso_model,type="coefficients")
lasso.coef[lasso.coef!=0]
lasso_nonzerocoef <- predict(lasso_model,type="nonzero")
colnames(data_qsar[,lasso_nonzerocoef$s0+1])
# model evaluation on lasso model using all non-zero descriptors
lasso_fittings <- predict(lasso_model, s = lambda_best_lasso, newx = x_train)
lasso_predictions <- predict(lasso_model, s = lambda_best_lasso, newx =
x_test)
evaluate_results(y_test, lasso_predictions, test_set)
evaluate_results(y_train, lasso_fittings, training_set)
```

    [1] 100



    
![png](Feature_selection_files/Feature_selection_11_1.png)
    



<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>-16.0130462430869</li><li>-0.00352928146601359</li><li>0.148827432572551</li><li>-0.628292587502564</li><li>-3.84326111653081e-07</li><li>-0.0743499949304882</li><li>0.291162379572836</li><li>-0.0437621280442267</li><li>-3.82883797083726e-06</li><li>-0.00396459114774324</li><li>-138.867203265117</li><li>2.7949620538661</li><li>-9.71674425865623</li><li>-4.16160352865241</li><li>-0.493074218693939</li><li>-0.197070283389385</li><li>-0.148983099209511</li><li>-0.00603386041367994</li><li>0.0015753452657488</li><li>-2.09494358340313e-05</li><li>19.3345183000582</li><li>-0.00742736451032618</li><li>-0.000431155658331097</li><li>-0.0313200631383102</li><li>0.00292892129242884</li><li>-0.00951027529792586</li><li>0.0418998156858545</li><li>0.13040122352195</li><li>-0.00779441737304536</li><li>-0.0914561200706204</li><li>-0.177942270528242</li><li>-0.002394365706595</li><li>0.225796648255782</li><li>0.539449844265398</li><li>0.840607514668115</li><li>0.0234395459941934</li></ol>




<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'AM1_HF'</li><li>'ast_violation'</li><li>'b_max1len'</li><li>'DCASA'</li><li>'dipoleY'</li><li>'dipoleZ'</li><li>'E_ang'</li><li>'E_rnb'</li><li>'E_rsol'</li><li>'E_strain'</li><li>'GCUT_PEOE_0'</li><li>'GCUT_PEOE_1'</li><li>'GCUT_PEOE_2'</li><li>'GCUT_SLOGP_0'</li><li>'h_pKb'</li><li>'h_pstates'</li><li>'PEOE_VSA.0'</li><li>'PEOE_VSA.1.1'</li><li>'PEOE_VSA_POS'</li><li>'petitjean'</li><li>'SlogP_VSA3'</li><li>'SlogP_VSA8'</li><li>'SMR_VSA4'</li><li>'SMR_VSA7'</li><li>'vsa_acc'</li><li>'vsa_other'</li><li>'vsurf_A'</li><li>'vsurf_DD13'</li><li>'vsurf_DD23'</li><li>'vsurf_DW12'</li><li>'vsurf_DW13'</li><li>'vsurf_ID1'</li><li>'vsurf_ID3'</li><li>'vsurf_IW2'</li><li>'vsurf_IW7'</li></ol>




<table class="dataframe">
<caption>A data.frame: 12 × 2</caption>
<thead>
	<tr><th></th><th scope=col>lnKD</th><th scope=col>lnKD.1</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>2</th><td>0.34456179</td><td>-2.131984e+00</td></tr>
	<tr><th scope=row>6</th><td>0.22400086</td><td>-4.747976e-01</td></tr>
	<tr><th scope=row>10</th><td>0.01475028</td><td> 9.998310e-01</td></tr>
	<tr><th scope=row>13</th><td>0.04390299</td><td> 9.913746e-01</td></tr>
	<tr><th scope=row>19</th><td>0.25568366</td><td> 5.427119e-01</td></tr>
	<tr><th scope=row>24</th><td>0.48732364</td><td> 8.598265e-01</td></tr>
	<tr><th scope=row>27</th><td>0.20593396</td><td>-2.150123e+04</td></tr>
	<tr><th scope=row>32</th><td>0.96412899</td><td> 4.833237e-01</td></tr>
	<tr><th scope=row>34</th><td>0.31965730</td><td>-3.522995e+01</td></tr>
	<tr><th scope=row>36</th><td>0.41677210</td><td> 3.781934e-01</td></tr>
	<tr><th scope=row>39</th><td>0.04886899</td><td> 9.917872e-01</td></tr>
	<tr><th scope=row>48</th><td>0.96748482</td><td>-1.337226e+00</td></tr>
</tbody>
</table>




<table class="dataframe">
<caption>A data.frame: 36 × 2</caption>
<thead>
	<tr><th></th><th scope=col>lnKD</th><th scope=col>lnKD.1</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>1</th><td>0.025316068</td><td>0.9988648</td></tr>
	<tr><th scope=row>3</th><td>0.014286497</td><td>0.9845556</td></tr>
	<tr><th scope=row>4</th><td>0.001797751</td><td>0.9990036</td></tr>
	<tr><th scope=row>5</th><td>0.025991215</td><td>0.9607456</td></tr>
	<tr><th scope=row>7</th><td>0.010787365</td><td>0.9965625</td></tr>
	<tr><th scope=row>8</th><td>0.025305323</td><td>0.8781310</td></tr>
	<tr><th scope=row>9</th><td>0.007877363</td><td>0.9954080</td></tr>
	<tr><th scope=row>11</th><td>0.022760976</td><td>0.8130770</td></tr>
	<tr><th scope=row>12</th><td>0.016993506</td><td>0.9995385</td></tr>
	<tr><th scope=row>14</th><td>0.019230594</td><td>0.9988566</td></tr>
	<tr><th scope=row>15</th><td>0.033135378</td><td>0.9987086</td></tr>
	<tr><th scope=row>16</th><td>0.009205695</td><td>0.9974306</td></tr>
	<tr><th scope=row>17</th><td>0.006232780</td><td>0.9993132</td></tr>
	<tr><th scope=row>18</th><td>0.016281516</td><td>0.9953077</td></tr>
	<tr><th scope=row>20</th><td>0.023231315</td><td>0.9136957</td></tr>
	<tr><th scope=row>21</th><td>0.019958698</td><td>0.9976389</td></tr>
	<tr><th scope=row>22</th><td>0.002884242</td><td>0.9999717</td></tr>
	<tr><th scope=row>23</th><td>0.008923600</td><td>0.9998404</td></tr>
	<tr><th scope=row>25</th><td>0.018165764</td><td>0.9907126</td></tr>
	<tr><th scope=row>26</th><td>0.005756368</td><td>0.9996894</td></tr>
	<tr><th scope=row>28</th><td>0.005015090</td><td>0.9997078</td></tr>
	<tr><th scope=row>29</th><td>0.006542657</td><td>0.9997362</td></tr>
	<tr><th scope=row>30</th><td>0.028536696</td><td>0.9968921</td></tr>
	<tr><th scope=row>31</th><td>0.007806772</td><td>0.9989256</td></tr>
	<tr><th scope=row>33</th><td>0.033356390</td><td>0.9978950</td></tr>
	<tr><th scope=row>35</th><td>0.017113126</td><td>0.9795642</td></tr>
	<tr><th scope=row>37</th><td>0.005802401</td><td>0.9998232</td></tr>
	<tr><th scope=row>38</th><td>0.049930739</td><td>0.9774642</td></tr>
	<tr><th scope=row>40</th><td>0.003203095</td><td>0.9963653</td></tr>
	<tr><th scope=row>41</th><td>0.015493340</td><td>0.9832864</td></tr>
	<tr><th scope=row>42</th><td>0.009914384</td><td>0.9930038</td></tr>
	<tr><th scope=row>43</th><td>0.017011834</td><td>0.9701513</td></tr>
	<tr><th scope=row>44</th><td>0.016359827</td><td>0.9984996</td></tr>
	<tr><th scope=row>45</th><td>0.004229424</td><td>0.9998920</td></tr>
	<tr><th scope=row>46</th><td>0.028836635</td><td>0.9813104</td></tr>
	<tr><th scope=row>47</th><td>0.001100779</td><td>0.9999914</td></tr>
</tbody>
</table>




    
![png](Feature_selection_files/Feature_selection_11_6.png)
    



```R
print (append(lasso_nonzerocoef$s0+1, 1, 0))

data_step <- training_set[,append(lasso_nonzerocoef$s0+1, 1, 0)]
```

     [1]   1   6  13  33  39  44  45  46  51  52  56  64  65  66  67  80  81 100 108
    [20] 117 118 135 138 143 145 149 152 153 161 162 163 164 173 174 178 183



```R
m <- 3
idx <- combn(rep(1:(length(data_step)-1)),m)

```


```R
# Your R code here
suppressWarnings({
results <- NULL
for (i in 1:ncol(idx)) {
data_exhau <- data_step[,append(idx[,i]+1,1,0)] 
mdl_exhau <- lm(lnKD~.,data=data_exhau)
predict <- predict(mdl_exhau,newdata = test_set)
fitted <- mdl_exhau$fitted.values
a <- evaluate_results(test_set$lnKD,predict,test_set)
b <- evaluate_results(training_set$lnKD,fitted,training_set)

    result <- data.frame(test=a, 
                     train=b
) 
    results <- rbind(results,result)
}
# idrows find all candidates with top performance, and print out the model summary for statistical significance check
})

```


```R

#idrows <- which(results$test.R_squared>=0.70&results$train.R_squared>=.70)

#for (val in idrows) {
#data_exhau <- data_step[,append(idx[,val-1]+1,1,0)] 
#    mdl_exhau <- lm(lnKD~.,data=data_exhau)
    
#s <- summary(mdl_exhau)
# print(s)
# print(val)
# cat("R2_test:", results[val,2])
#}

#s <- summary(mdl_exhau) 
#print(s)
#print(val)
#print("R2_test:", results[val,2])

# plot the curve for the top model
library(ggplot2)
# load the model
mdl <- lm(formula = "lnKD~1+PEOE_VSA_POS+vsurf_DW12+vsa_other+vsurf_ID3", data = training_set) 
summary(mdl)
predict <- predict(mdl,newdata = data_qsar) 
id <- numeric(48)
id[-trainid] <- 1
data_plot <- cbind(predict,data_qsar$lnKD,id)
colnames(data_plot) <- c("predict", "obs","id")
ggplot(as.data.frame(data_plot), aes(x=obs,y=predict))+ ggtitle(expression("Baseline model of lnK"[D]*"")) + xlab(expression("Observed lnK"[D]*"")) + ylab(expression("Predicted
lnK"[D]*""))+
# THE DATA POINT
geom_point(aes(color = factor(id)),size = 5,alpha =1) + xlim(min(data_qsar$lnKD)-2,max(data_qsar$lnKD)+2)+ ylim(min(data_qsar$lnKD)-2,max(data_qsar$lnKD)+2)+
scale_color_manual(labels = c("Training set", "Test set"), values =
c("dodgerblue", "red2"))+
# title
theme_bw()+ theme(axis.ticks.length=unit(.4,"lines"))+ theme(panel.grid.major = element_blank(),
panel.grid.minor = element_blank())+ theme(axis.text.y = element_text(size = 20),
axis.text.x = element_text(size=20),
axis.title = element_text(size = 25,face = 'bold'),title =element_text(size = 25,face = 'bold') )+
# legend
theme(legend.title = element_blank())+
theme(legend.text = element_text(colour="black", size=20, face="bold"))+ theme(legend.position = c(0.80, 0.1))+
# rec
theme(panel.background = element_rect(colour = "black", size = 3.5))+ # ref line
geom_abline(intercept = 0, slope = 1, color="black",
linetype="dashed", size=1.5) 
ggsave("KDmdl.tiff", units="in", width=8, height=8, dpi=600)
    



```


    
    Call:
    lm(formula = "lnKD~1+PEOE_VSA_POS+vsurf_DW12+vsa_other+vsurf_ID3", 
        data = training_set)
    
    Residuals:
        Min      1Q  Median      3Q     Max 
    -2.7574 -0.6231  0.1389  0.5923  2.4163 
    
    Coefficients:
                  Estimate Std. Error t value Pr(>|t|)    
    (Intercept)  -10.03906    0.88572 -11.334 1.48e-12 ***
    PEOE_VSA_POS  -0.01521    0.00232  -6.557 2.54e-07 ***
    vsurf_DW12    -0.37242    0.12798  -2.910  0.00663 ** 
    vsa_other      0.05441    0.01029   5.286 9.46e-06 ***
    vsurf_ID3      1.68556    0.46762   3.605  0.00108 ** 
    ---
    Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
    
    Residual standard error: 1.237 on 31 degrees of freedom
    Multiple R-squared:  0.7686,	Adjusted R-squared:  0.7387 
    F-statistic: 25.74 on 4 and 31 DF,  p-value: 1.817e-09



    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “font metrics unknown for character 0xa”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “font metrics unknown for character 0xa”
    Warning message in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    “font metrics unknown for character 0xa”
    Warning message in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    “font metrics unknown for character 0xa”



    
![png](Feature_selection_files/Feature_selection_15_2.png)
    



```R

```
