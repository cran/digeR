\name{digeR}
\Rdversion{1.1}
\alias{digeR}
\title{digeR GUI 
}
\description{ Start the Graphical User Interface for digeR. \cr
digeR supports spots correlation analysis, score plot, classification, feature selection and power analysis.
}
\usage{
digeR()
}
\details{
\strong{digeR GUI options:}
      
\tabular{lll}{
    \tab \strong{File}                       \tab Read in data and image, quit \cr
    \tab Open	                               \tab upload the txt file \cr
    \tab Upload_gel_image	                   \tab upload the JPG image as a reference for spots correlation analysis \cr
    \tab Quit	                               \tab dispose the GUI \cr
    \tab --------------------\tab \cr              
    \tab \strong{Correlation}               \tab  Spots correlation analysis    \cr
    \tab Dataset					                  \tab  select the group to look at   \cr
    \tab Spot List				                  \tab  select the spot to look at   \cr
    \tab Selected feature		                \tab	upload the feature list from feature selection  \cr
    \tab Load features				              \tab  upload the feature list from an saved R workspace \cr
    \tab Show the correlation	              \tab  plot the spots with required correlation \cr
    \tab Correlation Coefficiency	          \tab	change the coefficiency threshold \cr
    \tab Show spot ID			                  \tab	plot spots with ID \cr
    \tab Show number 			                  \tab  Show ID for those spots with required correlation \cr
    \tab --------------------\tab  \cr              
    \tab \strong{Score Plot}                         \tab PCA and PLSR score plot \cr
    \tab Plot Type						             \tab select either PCA or PLSR score plot \cr
    \tab Top N component				           \tab plot score plot with top N components \cr
    \tab Pairwise							           \tab plot selected 2 components \cr
    \tab Component 1 and 2				        \tab two components in the pairwise plot \cr
    \tab Group							              \tab set the color for the two groups \cr
    \tab With label						          \tab plot the sample ID  \cr
    \tab Scaling							            \tab scale the data before plotting \cr
    \tab Plot								            \tab plot the score plot  \cr
    \tab --------------------\tab \cr             
    \tab \strong{Classification}                  \tab Classification     \cr
    \tab Methods						              \tab select the method for the classification \cr
    \tab Scaling							            \tab scale the data before classification \cr
    \tab \emph{Arguments} \tab                      \cr
    \tab Method							            \tab way for estimate the covariance matrix. \cr
    \tab "moment" 						            \tab standard estimators of the mean and variance \cr
    \tab "mle"							              \tab MLEs,  \cr
    \tab "mve" 							             \tab to use cov.mve \cr
    \tab "t" 		 						            \tab robust estimates based on a t distribution \cr
    \tab nComp							              \tab number of component for fitting PCR or PLSR \cr
    \tab N fold CV						            \tab number of fold in the cross validation  \cr
    \tab nboot							              \tab number of bootstrap in the classification  \cr
    \tab Selected feature					       \tab upload the feature list from feature selection \cr
    \tab Load features					          \tab upload the feature list from an saved R workspace \cr
    \tab leave one out cv					        \tab classification with leave one out cross validation \cr
    \tab N fold cv						            \tab classification with n fold cross validation \cr
    \tab Bootstrap						            \tab classification with bootstrap \cr
    \tab Run classification				      \tab press button to do the classification  \cr
    \tab Save							              \tab save the prediction results into an R workspace  \cr
    \tab Legened							            \tab where the legend will be put \cr
    \tab ROC curve						            \tab generate ROC plot          \cr
    \tab Prediction result					      \tab store the classification results in the selected items \cr
    \tab --------------------\tab \cr               
    \tab \strong{Feature Selection}               \tab Select important features \cr
    \tab Method							            \tab select feature selection method  \cr
    \tab Scaling							            \tab scale the data before feature selection  \cr
    \tab \emph{Arguments} \tab                      \cr
    \tab Method							            \tab same as Method in Classification \cr
    \tab Ncomp							              \tab same as ncomp in Classification  \cr
    \tab Top								              \tab select the top n variables from the feature selection \cr
    \tab Ntree							              \tab number of trees to grow in randomForest \cr
    \tab Mfinal	                        \tab the number of iterations for which boosting is run or the number of trees to use \cr
    \tab Run feature selection				    \tab press to start feature selection  \cr
    \tab Select featuers					        \tab store the selected features in the selected items \cr
    \tab Save features					          \tab save the features into an R workspace  \cr
    \tab --------------------\tab \cr                 
    \tab \strong{Power}                   \tab Power analysis      \cr
    \tab Single Spots						         \tab univariate power analysis   \cr
    \tab Gel								              \tab multivariate power analysis for experiment design  \cr
    \tab Significant level				        \tab set the significant level      \cr
    \tab Power							              \tab set the power level to be achieved  \cr
    \tab Sample size per group	          \tab sample size for achieving certain significant level and power in each group \cr
    \tab Spot Number					            \tab set the spots to calculated \cr
    \tab Calculate 	                    \tab calculate the one being left blank (either power, sample size or significant level) \cr
}

\strong{The details can be found in digeR manual}
}

\author{
Yue Fan \email{yue.fan@ucd.ie},Thomas Brendan Murphy, R. William G. Watson
}
\note{
digeR is built upon gWidgets package. Make sure gWidgets package is properly installed.

}


