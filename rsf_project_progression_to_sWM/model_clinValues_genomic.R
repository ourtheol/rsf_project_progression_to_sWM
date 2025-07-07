library(randomForestSRC)
library(pROC)
library(survival)
library(ggsurvfit)
library(dplyr)

setwd("./")
 
# there is randomization!
set.seed(1)

d <- read.table("./data/training_cohort_initial64patients.tsv", header = TRUE, row.names = 1)
print(head(d))


# # Use our RSFclin+genomic model or uncomment the following and create your own model
obj <- readRDS("./models/obj_model3_clinValues_genomic.RDS")
obj 



# # Find the optimal mtry and nodesize tuning parameter for a random forest using out-of-bag (OOB) error.
# t <- tune(Surv(Months,Progression) ~ ., d, ntree = 1000, nsplit = 50)
# t$optimal
# # The following features were used when constructing our RSFclin+genomic (obj_model3_clinValues_genomic.RDS)
# # nodesize     mtry
# # 10           24
# 
# 
# # number of trees
# ntree = 1000
# # size of terminal node
# nodesize = t$optimal[[1]]
# # number of variables randomly selected as candidates for splitting a node
# mtry = t$optimal[[2]]
# # fixed number of randomly selected split-points, to reduce computations when considering all possible split-values for a variable,
# # this also mitigates the well known tree bias of favoring splits on variables with a large number of split-points (continuous variables) or factors with a large number of categorical labels
# nsplit = 50
# 
# 
# # Plot cumulative OOB error rates as a function of number of trees and variable importance (VIMP).
# # To get OOB error rates for every tree, use the option block.size = 1 when growing or restoring the forest.
# # Likewise, to view VIMP, use the option importance when growing or restoring the forest.
# pdf("./results/cumulative_error_VIMP.pdf", 30, 30)
# plot(rfsrc(Surv(Months, Progression) ~ ., d, block.size = 1, ntree = ntree, nodesize = nodesize, mtry = mtry, nsplit = nsplit, importance = TRUE))
# dev.off()
# 
# 
# obj <- rfsrc(Surv(Months,Progression)~., d, ntree = ntree, nodesize = nodesize, mtry = mtry, nsplit = nsplit,  importance = TRUE)
# 
# # Save the model
# saveRDS(obj, "./models/obj_model3_clinValues_genomic_new.RDS")



# Get the risk scores for each patient
scores.table <- as.data.frame(cbind(rownames(d),obj$predicted.oob))
colnames(scores.table) <- c("Sample", "score")
scores.table$score <- as.numeric(scores.table$score)

write.table(scores.table,"./results/risk_scores.csv", quote = FALSE, row.names = FALSE)



# Info df for each patient
patient.states <- read.table("./data/patient_states_training_cohort.tsv", header = TRUE)

sc <- merge(scores.table, patient.states, by = "Sample")
head(sc)
str(sc)

# Calculate the distance matrix
dist_matrix <- dist(sc$score, method = "euclidean")

# Perform hierarchical clustering
hclust_result <- hclust(dist_matrix, method = "ward.D")

# Plot the dendrogram with patient names
pdf("./results/hist_scores.pdf")
plot(hclust_result, main = "Hierarchical Clustering Dendrogram",
     xlab = "Patients", sub = "", labels = sc$Sample, cex = 0.8)
dev.off()

# Plot a scatter plot of the scores per patient (y-axis) alonng the months they have been observed
# dots are colored according to the clinical state per patient
scatter <- ggplot(sc, aes(x=Months, y=score, color=State, label=Sample)) + 
                  geom_point(size=2) +
                  scale_color_manual(values=c( 'gold',"deepskyblue3",  'forestgreen' )) +
                  geom_point() +
                  geom_text() +
                  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1, size=18), axis.text.y = element_text(size=16),
                        legend.key.size = unit(2, 'cm'),legend.text = element_text(size=15),legend.title = element_text(size=0),
                        axis.title=element_text(size=17,face="bold"),
                        panel.background = element_rect(fill = "white", colour = "black",linewidth = 1.5, linetype = "solid"),
                        panel.grid.major = element_line(linewidth = 0.5, linetype = 'solid', colour = "grey"),
                        aspect.ratio=1,
                        legend.position = "none") 

pdf("./results/scatter_scores.pdf")
print(scatter)
dev.off()



# ROC-analysis for patient stratification into groups according to RSF scores
roc_obj <- roc(sc$Progression, sc$score, direction = "auto")  # "<" means low score = high risk
plot(roc_obj, main = paste("AUC =", round(auc(roc_obj), 3)))

# Best threshold based on Youden’s index (maximizes Sensitivity + Specificity - 1)
best_thresh <- coords(roc_obj, "best", ret = "threshold")
message(paste0("ROC-analysis for patient stratification threshold:", best_thresh))
# the threshold in our RSFclin+genomic model is 1.09

# Assign patients into groups according to the RSF scrores
sc <- sc %>%
  mutate(score_group = ifelse(score >= best_thresh[[1]], "score_group_1", "score_group_2"))

write.table(sc, "./results/score_groups.csv", row.names = FALSE, quote = FALSE)



jk.obj <- subsample(obj)
pdf("./results/VIMPsur.pdf", width = 25, height = 30)
par(oma = c(0.5, 10, 0.5, 0.5))
par(cex.axis = 2.0, cex.lab = 2.0, cex.main = 2.0, mar = c(6.0,17,1,1), mgp = c(4, 1, 0))
plot(jk.obj, xlab = "Variable Importance (x 100)", cex = 1.2)
dev.off()

# find interactions
intr <- find.interaction(obj, method = "vimp", nvar = 20)
write.table(intr,"./results/interactions.csv", sep = ",", quote = FALSE, row.names = FALSE)

# Extract maximal subtree information from a RF-SRC object. Used for variable selection and identifying interactions between variables.
v.max <- max.subtree(obj)

# first and second order depths
print(round(v.max$order, 3))

# the minimal depth is the first order depth
print(round(v.max$order[, 1], 3))

# strong variables have minimal depth less than or equal to the following threshold
print(v.max$threshold)

# this corresponds to the set of variables
print(v.max$topvars)

write.table(print(v.max$topvars),"./results/topvars.csv", quote = FALSE, row.names = FALSE)



# Minimal depth variable selection
# use larger node size which is better for minimal depth

# default call corresponds to minimal depth selection
vs.pbc <- var.select(object = obj)
topvars <- vs.pbc$topvars

# the above is equivalent to
max.subtree(obj)$topvars

# different levels of conservativeness
var.select(object = obj, conservative = "low")
var.select(object = obj, conservative = "medium")
var.select(object = obj, conservative = "high")

write.table(topvars,"./results/topvars2.csv", quote = FALSE, row.names = FALSE)



# #extract a single tree from a forest and plot it on your browser
# plot(get.tree(obj, 5)) #library(data.tree) is imported



# Use a subset to focus on specific individuals
pdf("./results/surv.pdf")
plot.survival(obj, subset = 1) # sample T_1546
# all progressed together
plot.survival(obj, subset = c(1,55,56,57,59,61,62,63,64))
dev.off()

# Select subset of patients and make custom cumulative incidence plots based on the OOB predictions
# This works with the OOB predictions
subset = c(1,55,56,57,59,61,62,63,64)
title = "aWM_pr"

brier.obj <- get.brier.survival(obj, subset = subset, cens.model=c("km", "rfsrc"))  
surv.ensb <- brier.obj$surv.ensb
event.info <- brier.obj$event.info
surv.mean.ensb <- rowMeans(surv.ensb, na.rm = TRUE)

pdf("results/predicted_CI_aWMpr.pdf")
matplot(obj$time.interest,
        surv.ensb,
        xlab = "Time (Months)",
        ylab = "OOB Survival",
        main = title,
        type = "l",
        col = 1,
        lty = 3) 
lines(event.info$time.interest, surv.mean.ensb, lty = 1, col = 2, lwd = 3)
dev.off()



# Validation
val <- read.table("./data/validation_cohort.tsv", header = TRUE, row.names = 1)

y.pred <- predict(obj, newdata = val)
# predicted scores
y.pred$predicted

scores.table.validation <- as.data.frame(cbind(rownames(val),y.pred$predicted))
colnames(scores.table.validation) <- c("Sample", "score")
scores.table.validation$score <- as.numeric(scores.table.validation$score)



# Assign patients into groups according to the RSF scrores
scores.table.validation <- scores.table.validation %>%
  mutate(score_group = ifelse(score >= best_thresh[[1]], "score_group_1", "score_group_2"))
          
patient.states.val <- read.table("./data/patient_states_validation_cohort.tsv", header = TRUE)

scores.table.validation <- merge(scores.table.validation, patient.states.val, by = "Sample")

write.table(scores.table.validation, "./results/score_groups_validation.csv", quote = FALSE, row.names = FALSE)



# Survival analysis for the training cohort, sc table
# Conduct between-group significance tests using a log-rank test, to compare survival times between groups.
formula <- Surv(Months, Progression) ~ score_group

fit <- survfit(formula = formula, data = sc)

diff <- survdiff(formula = formula, data = sc)

# Cox proportional hazards model to illustrate the estimated probability of survival over time.
fit.cox <- coxph(formula = formula, data = sc)

c <- concordance(fit.cox, timewt="n/G2") # Uno's weighting


# Write all outputs to the same file
output_file <- "./results/survival_training_cohort.txt"
sink(output_file)

cat("==== Formula ====\n")
print(formula)

cat("==== Kaplan-Meier Fit (survfit) ====\n")
print(fit)

cat("\n==== Log-rank Test (survdiff) ====\n")
print(diff)

cat("\n==== Cox Proportional Hazards Model ====\n")
print(summary(fit.cox))  # use summary for more detailed output

cat("\n==== Concordance Index ====\n")
print(c)

# Stop redirecting output
sink()


p <- survfit2(formula = formula, data = sc) %>%  
  ggsurvfit(type = "risk") +  #type = "risk" when you want to plot the cumulative
  labs(
    x = "Months",
    y = "Cumulative incidence"   #"Overall survival probability" for the survival curves
  ) +
  add_confidence_interval() +
  add_risktable() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1, size=18), axis.text.y = element_text(size=16),
        legend.key.size = unit(2, 'cm'),legend.text = element_text(size=15),legend.title = element_text(size=0),
        axis.title=element_text(size=17,face="bold"),
        panel.background = element_rect(fill = "white", colour = "black",linewidth = 1.5, linetype = "solid"),
        panel.grid.major = element_line(linewidth = 0.5, linetype = 'solid', colour = "grey"),
        aspect.ratio=1, legend.position = "none")



# Survival analysis for the validation cohort, scores.table.validation table.
# Conduct between-group significance tests using a log-rank test, to compare survival times between groups.
fit <- survfit(formula = formula, data = scores.table.validation)

diff <- survdiff(formula = formula, data = scores.table.validation)

# Cox proportional hazards model to illustrate the estimated probability of survival over time.
fit.cox <- coxph(formula = formula, data = scores.table.validation)

c <- concordance(fit.cox, timewt="n/G2") # Uno's weighting


# Write all outputs to the same file
output_file <- "./results/survival_validation_cohort.txt"
sink(output_file)

cat("==== Formula ====\n")
print(formula)

cat("==== Kaplan-Meier Fit (survfit) ====\n")
print(fit)

cat("\n==== Log-rank Test (survdiff) ====\n")
print(diff)

cat("\n==== Cox Proportional Hazards Model ====\n")
print(summary(fit.cox))  # use summary for more detailed output

cat("\n==== Concordance Index ====\n")
print(c)

# Stop redirecting output
sink()


# Survival curves
p2 <- survfit2(formula = formula, data = scores.table.validation) %>%  
  ggsurvfit(type = "risk") +  #type = "risk" when you want to plot the cumulative
  labs(
    x = "Months",
    y = "Cumulative incidence"   #"Overall survival probability" for the survival curves
  ) +
  add_confidence_interval() +
  add_risktable() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1, size=18), axis.text.y = element_text(size=16),
        legend.key.size = unit(2, 'cm'),legend.text = element_text(size=15),legend.title = element_text(size=0),
        axis.title=element_text(size=17,face="bold"),
        panel.background = element_rect(fill = "white", colour = "black",linewidth = 1.5, linetype = "solid"),
        panel.grid.major = element_line(linewidth = 0.5, linetype = 'solid', colour = "grey"),
        aspect.ratio=1, legend.position = "none")


pdf("./results/survival_rsf_score_groups.pdf")
p
p2
dev.off()