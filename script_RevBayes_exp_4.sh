#!/bin/bash

calls="alarm"

# Generating the Q matrix using R
Rscript script_write_matrix_sig_exp_4.R $calls 

# Fomating the matrix file
sed -ie 's/$/],/g' Q_tmp
sed -ie 's/^/[/g' Q_tmp
sed -ie '1s/^/[/' Q_tmp
sed -ie '$s/\(.*\),/\1]/g' Q_tmp
rm Q_tmpe

p=$( cat param_tmp )
np=$( cat n_param_tmp )
rm param_tmp
rm n_param_tmp

# Begining of the mcmc_primates.Rev new file
echo '# Automatically generated mcmc_primates file for RevBayes' > mcmc_primates_tmp.Rev

echo 'monitors = VectorMonitors()' >> mcmc_primates_tmp.Rev
echo 'moves    = VectorMoves()' >> mcmc_primates_tmp.Rev

echo "mu13 ~ dnExponential( 1.0 )" >> mcmc_primates_tmp.Rev
echo "mu13.setValue( 0.1 )" >> mcmc_primates_tmp.Rev

echo "mu99 ~ dnExponential( 1.0 )" >> mcmc_primates_tmp.Rev
echo "mu99.setValue( 0.1 )" >> mcmc_primates_tmp.Rev

echo "mu_0 ~ dnExponential( 1.0 )" >> mcmc_primates_tmp.Rev
echo "mu_0.setValue( 0.1 )" >> mcmc_primates_tmp.Rev

echo "lambda_0 ~ dnUniform(0, 10)" >> mcmc_primates_tmp.Rev
echo "lambda_0.setValue(1)" >> mcmc_primates_tmp.Rev

echo "delay_0 ~ dnUniform(0, $np + 2)" >> mcmc_primates_tmp.Rev
echo "delay_0.setValue( 3 )" >> mcmc_primates_tmp.Rev

for i in $p
do
	#echo "mu$i := mu_0 * exp( - ($i-1) / (25*lambda_0) )" >> mcmc_primates_tmp.Rev
	echo "mu$i := mu_0 / (1 + exp(lambda_0 * ($i - delay_0) ) )" >> mcmc_primates_tmp.Rev
done


echo 'Q_morpho :=' >> mcmc_primates_tmp.Rev
cat Q_tmp >> mcmc_primates_tmp.Rev

rm Q_tmp

echo "moves.append( mvScale(mu_0,lambda=1, weight=1.0) )" >> mcmc_primates_tmp.Rev
echo "moves.append( mvScale(mu13,lambda=1, weight=1.0) )" >> mcmc_primates_tmp.Rev
echo "moves.append( mvScale(mu99,lambda=1, weight=1.0) )" >> mcmc_primates_tmp.Rev
echo "moves.append( mvScale(lambda_0,lambda=1, weight=1.0) )" >> mcmc_primates_tmp.Rev
echo "moves.append( mvScale(delay_0,lambda=1, weight=1.0) )" >> mcmc_primates_tmp.Rev

echo "morpho <- readDiscreteCharacterData(\"lexicons_primates_$calls.nex\")"  >> mcmc_primates_tmp.Rev
echo "Q := fnFreeK(Q_morpho, FALSE)" >> mcmc_primates_tmp.Rev
echo "phylogeny <- readTrees(\"primatesx100.tre\")[1]" >> mcmc_primates_tmp.Rev
echo "phyMorpho ~ dnPhyloCTMC(tree=phylogeny, Q=Q, type=\"Standard\")" >> mcmc_primates_tmp.Rev
echo "phyMorpho.clamp(morpho)" >> mcmc_primates_tmp.Rev
echo "mymodel = model(phylogeny)" >> mcmc_primates_tmp.Rev
echo "monitors.append( mnModel(filename=\"output_sig_cat_$calls/mk_$calls.log\", printgen=1) )" >> mcmc_primates_tmp.Rev
echo "monitors.append( mnScreen(printgen=100) )" >> mcmc_primates_tmp.Rev
echo "mymcmc = mcmc(mymodel, monitors, moves, nruns=1, combine=\"mixed\")" >> mcmc_primates_tmp.Rev
echo "mymcmc.run(generations=10000, tuningInterval=1000)" >> mcmc_primates_tmp.Rev
echo "mymcmc.operatorSummary()" >> mcmc_primates_tmp.Rev
echo "q()" >> mcmc_primates_tmp.Rev

# Running RavBayes
./rb mcmc_primates_tmp.Rev

rm mcmc_primates_tmp.Rev