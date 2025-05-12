library(randomForestSRC)
library(pROC)


setwd("./")
 
# there is randomization!
set.seed(1)

d <- read.table("./data/training_cohort_initial64patients.tsv", header = TRUE, row.names = 1)
print(head(d))


# Find the optimal mtry and nodesize tuning parameter for a random forest using out-of-bag (OOB) error.
t <- tune(Surv(Months,Progression) ~ ., d, ntree = 1000, nsplit = 50)
t$optimal
# The following features were used when constructing our RSFclin+genomic (obj_model3_clinValues_genomic.RDS)
# nodesize     mtry  
# 10           24


# number of trees
ntree = 1000
# size of terminal node
nodesize = t$optimal[[1]]
# number of variables randomly selected as candidates for splitting a node
mtry = t$optimal[[2]]
# fixed number of randomly selected split-points, to reduce computations when considering all possible split-values for a variable, 
# this also mitigates the well known tree bias of favoring splits on variables with a large number of split-points (continuous variables) or factors with a large number of categorical labels
nsplit = 50


# Plot cumulative OOB error rates as a function of number of trees and variable importance (VIMP).
# To get OOB error rates for every tree, use the option block.size = 1 when growing or restoring the forest. 
# Likewise, to view VIMP, use the option importance when growing or restoring the forest.
pdf("./results/cumulative_error_VIMP.pdf", 30, 30)
plot(rfsrc(Surv(Months, Progression) ~ ., d, block.size = 1, ntree = ntree, nodesize = nodesize, mtry = mtry, nsplit = nsplit, importance = TRUE))
dev.off()


obj <- rfsrc(Surv(Months,Progression)~., d, ntree = ntree, nodesize = nodesize, mtry = mtry, nsplit = nsplit,  importance = TRUE) 

# Save the model
saveRDS(obj, "./models/obj_model3_clinValues_genomic_new.RDS")


# Get the risk scores for each patient
scores.table <- as.data.frame(cbind(rownames(d),obj$predicted.oob))
colnames(scores.table) <- c("Sample", "score")

write.table(scores.table,"./results/risk_scores.csv")



# Info df for each patient
patient.states <- read.table("./data/patient_states_training_cohort.tsv")

sc <- merge(scores.table, patient.states, by = "Sample")


# Calculate the distance matrix
dist_matrix <- dist(sc$score, method = "euclidean")

# Perform hierarchical clustering
hclust_result <- hclust(dist_matrix, method = "ward.D")

# Plot the dendrogram with patient names
hist <- plot(hclust_result, main = "Hierarchical Clustering Dendrogram",
             xlab = "Patients", sub = "", labels = sc$Sample, cex = 0.8)

pdf("./results/hist_scores.pdf", 30, 30)
print(hist)
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

pdf("./results/scatter_scores.pdf", 30, 30)
print(scatter)
dev.off()


# ROC-analysis for patient stratification into groups according to RSF scores

roc_obj <- roc(sc$Progression, sc$score, direction = "auto")  # "<" means low score = high risk
plot(roc_obj, main = paste("AUC =", round(auc(roc_obj), 3)))

# Best threshold based on Youden’s index (maximizes Sensitivity + Specificity - 1)
best_thresh <- coords(roc_obj, "best", ret = "threshold")
# the threshold in our RSFclin+genomic model is 1.09

# Assign patients into groups according to the RSF scrores
sc <- sc %>%
  mutate(score_group = ifelse(score >= best_thresh[[1]], "score_group_1", "score_group_2"))

write.table(sc, "./results/score_groups.csv")


jk.obj <- subsample(obj)
pdf("./results/VIMPsur.pdf", width = 25, height = 30)
par(oma = c(0.5, 10, 0.5, 0.5))
par(cex.axis = 2.0, cex.lab = 2.0, cex.main = 2.0, mar = c(6.0,17,1,1), mgp = c(4, 1, 0))
plot(jk.obj, xlab = "Variable Importance (x 100)", cex = 1.2)
dev.off()

# find interactions
intr <- find.interaction(obj, method = "vimp", nvar = 20)
write.table(intr,"./results/interactions.csv")


# Extract maximal subtree information from a RF-SRC object. Used for variable selection and identifying interactions between variables.
v.max <- max.subtree(obj)

# first and second order depths
print(round(v.max$order, 3))

# the minimal depth is the first order depth
print(round(v.max$order[, 1], 3))

# strong variables have minimal depth less than or equal
# to the following threshold
print(v.max$threshold)

# this corresponds to the set of variables
print(v.max$topvars)

write.table(print(v.max$topvars),"./results/topvars.csv")



## Minimal depth variable selection
## use larger node size which is better for minimal depth

# default call corresponds to minimal depth selection
vs.pbc <- var.select(object = obj)
topvars <- vs.pbc$topvars

# the above is equivalent to
max.subtree(obj)$topvars

# different levels of conservativeness
var.select(object = obj, conservative = "low")
var.select(object = obj, conservative = "medium")
var.select(object = obj, conservative = "high")

write.table(topvars,"./results/topvars2.csv")


# # extract a single tree from a forest and plot it on your browser
# plot(get.tree(obj, 5)) #library(data.tree) is imported


## use subset to focus on specific individuals
pdf("./results/surv.pdf")
plot.survival(obj, subset = 1) # sample T_1546
# all progressed together
plot.survival(obj, subset = c(1,55,56,57,59,61,62,63,64))
plot.survival(obj, subset = c(1,55,56,57,59,61,62,63,64), collapse = TRUE)
dev.off()


# Validation:
val <- read.table("./data/validation_cohort.tsv", header = TRUE, row.names = 1)

y.pred <- predict(obj, newdata = val)
# predicted scores
y.pred$predicted