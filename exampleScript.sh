#Example data set is in the directory "data"

######FILES##############
#Read counts for each site and sample. Rows are the sites and columns are the samples. First column is the name of the site.
data="data/signal.txt"
#Prior information of the proportion of each cell type and sample. Rows are samples and columns are cell types.
pPrior="data/pPrior.txt"
#Total read counts in each sample
scaleChip="data/sizes.txt"

#####HYPERPARAMETERS#####
#Hyperparameters w0, lower and upper bounds for cell type specific binding affinities and the number of optimization cycles.
w0="4";lower="0.01";upper="10000";maxIter="10"

#####OUTPUT#############
outputDIR="output";mkdir $outputDIR
outputP="$outputDIR/max_p.txt";outputSignal="$outputDIR/max_signal.txt"


#This script estimates the cell type specific signals in ChIP samples.
Rscript estimateChIP.R $data $w0 $pPrior $scaleChip $lower $upper $maxIter $outputP $outputSignal 


############################################################################################################
######FILES##############
#Read counts in input samples 1000 bases around each candidate site
#Rows are the sites and columns are the samples. First column is the name of the site.
data1000="data/input_1000.txt"
#Total read counts in each Input sample
scaleInput="data/sizes_Input.txt"

#####OUTPUT#############
output1000="$outputDIR/signal_input_1000.txt"

#This script estimates the cell type specific signals in input samples
Rscript estimateInput.R $data1000 $w0 $outputP $scaleInput $lower $upper $maxIter $output1000


############################################################################################################
######FILES##############
#Read counts in input samples 5000 bases around each candidate site
#Rows are the sites and columns are the samples. First column is the name of the site.
data5000="data/input_5000.txt"

#####OUTPUT#############
output1000="$outputDIR/signal_input_1000.txt"
output5000="$outputDIR/signal_input_5000.txt"

#This script estimates the cell type specific signals in input samples
Rscript  estimateInput.R $data5000 $w0 $outputP $scaleInput $lower $upper $maxIter $output5000


############################################################################################################
######FILES##############
#Read counts in input samples 10000 bases around each candidate site
#Rows are the sites and columns are the samples. First column is the name of the site.
data10000="data/input_10000.txt"

#####OUTPUT#############
output10000="$outputDIR/signal_input_10000.txt"

#This script estimates the cell type specific signals in input samples
Rscript estimateInput.R $data10000 $w0 $outputP $scaleInput $lower $upper $maxIter $output10000


############################################################################################################
######FILES##############
#Selected candidate binding sites, where columns are: chromosome, start, end and name of the site similar to the data matrix.
peaks="data/selected_peaks.txt"

#####OUTPUT#############
#Pvalues are written in the file called pvalues.txt
output="outputDIR/pvalues.txt"

#This script outputs pvalues of each peak.
Rscript getPvalues.R  $outputSignal $output1000 $output5000 $output10000 $peaks $data $scaleChip $scaleInput $output


