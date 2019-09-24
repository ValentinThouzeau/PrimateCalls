#!/bin/bash

let mod=97
calls="alarm"

# Generating the Q matrix using R
Rscript script_write_matrix_exp_1.R $calls $mod

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

for i in $p
do
	echo "mu$i ~ dnExponential( 2.0 )" >> mcmc_primates_tmp.Rev
    echo "mu$i.setValue( 0.1 )" >> mcmc_primates_tmp.Rev
done


echo 'Q_morpho :=' >> mcmc_primates_tmp.Rev
cat Q_tmp >> mcmc_primates_tmp.Rev

rm Q_tmp

for i in $p
do
	echo "moves.append( mvScale(mu$i,lambda=1, weight=1.0) )" >> mcmc_primates_tmp.Rev
done

echo "morpho <- readDiscreteCharacterData(\"lexicons_primates_$calls.nex\")"  >> mcmc_primates_tmp.Rev
echo "Q := fnFreeK(Q_morpho, FALSE)" >> mcmc_primates_tmp.Rev
echo "phylogeny <- readTrees(\"primatesx100.tre\")[1]" >> mcmc_primates_tmp.Rev
echo "phyMorpho ~ dnPhyloCTMC(tree=phylogeny, Q=Q, type=\"Standard\")" >> mcmc_primates_tmp.Rev
echo "phyMorpho.clamp(morpho)" >> mcmc_primates_tmp.Rev
echo "mymodel = model(phylogeny)" >> mcmc_primates_tmp.Rev
#echo "monitors.append( mnModel(filename=\"output$mod/mk.log\", printgen=1) )" >> mcmc_primates_tmp.Rev
#echo "monitors.append( mnScreen(printgen=100) )" >> mcmc_primates_tmp.Rev
#echo "mymcmc = mcmc(mymodel, monitors, moves, nruns=1, combine=\"mixed\")" >> mcmc_primates_tmp.Rev
#echo "mymcmc.run(generations=5000, tuningInterval=200)" >> mcmc_primates_tmp.Rev
#echo "mymcmc.operatorSummary()" >> mcmc_primates_tmp.Rev
echo "pow_p = powerPosterior(mymodel, moves, monitors, \"output_bd_$calls/mk_$calls.$mod.ss\", cats=49)" >> mcmc_primates_tmp.Rev
echo "pow_p.burnin(generations=10000,tuningInterval=1000)" >> mcmc_primates_tmp.Rev
echo "pow_p.run(generations=1000)" >> mcmc_primates_tmp.Rev
echo "ss = steppingStoneSampler(file=\"output_bd_$calls/mk_$calls.$mod.ss\", powerColumnName=\"power\", likelihoodColumnName=\"likelihood\")" >> mcmc_primates_tmp.Rev 
echo "ss_res = ss.marginal()" >> mcmc_primates_tmp.Rev
echo "write(ss_res, filename = \"ss_res_$calls\")" >> mcmc_primates_tmp.Rev

echo "q()" >> mcmc_primates_tmp.Rev

# Running RavBayes
./rb mcmc_primates_tmp.Rev

rm mcmc_primates_tmp.Rev

