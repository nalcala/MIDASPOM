######################################################################################
#### Script to run an example inference of extinction and colonization parameters ####
######################################################################################


#### run MIDASPOM ####
# Example parameter estimation in a linear habitat with 100m segments and mean dispersal of 400m
bin/MIDASPOM.out -m 400 -d 100 -i examples/input/occupancies.txt -o examples/output/posterior.txt

# Example parameter estimation under in situ die-off hypothesis, for surveys starting 10 years after the disturbance event
bin/MIDASPOM_dieoff.out -a 10 -e 0.71 -c 0.52 -m 400 -d 100 -s 151 -i examples/input/occupancies.txt -o examples/output/lh_dieoff.txt

# Example parameter estimation under habitat loss hypothesis, for surveys starting 10 years after the disturbance event
bin/MIDASPOM_loss.out -a 10 -e 0.71 -c 0.52 -m 400 -d 100 -s 151 -i examples/input/occupancies.txt -o examples/output/lh_loss.txt

# Example prediction of future probability of extinction without management
bin/MIDASPOM_future.out -a 50 -m 400 -d 100 -i examples/input/occupancies.txt -q examples/output/posterior.txt -o examples/output/pext_future_nomanagement.txt 

# Example prediction of future probability of extinction with an additional source population of size 1 at distance 500m
bin/MIDASPOM_future.out -a 50 -m 400 -d 100 -i examples/input/occupancies.txt -q examples/output/posterior.txt -o examples/output/pext_future_source.txt -S 1 -s 500


####  R scripts ####
# plot the posterior distributions of e and c
Rscript Rscript/plot_posterior.R  examples/output/posterior.txt examples/output/posterior.pdf 0 1 101

# test hypotheses and plot the posterior distributions of their parameters
Rscript Rscript/plot_posterior_dieoff.R  examples/output/lh_dieoff.txt examples/output/posterior_dieoff.pdf 0.1 100 151
Rscript Rscript/plot_posterior_loss.R  examples/output/lh_loss.txt examples/output/posterior_loss.pdf 0.1 100 200 4000 151 20

# perform hypothesis testing
Rscript Rscript/hypothesis_test.R  examples/output/lh_dieoff.txt examples/output/lh_loss.txt examples/output/hypotheses_test.pdf 8 0.1 100 151 20

# plot the probability of future extinction under different management scenarios
Rscript Rscript/plot_extinction.R  examples/output/pext_future_nomanagement.txt examples/output/pext_future_nomanagement.pdf 50 10000
Rscript Rscript/plot_extinction.R  examples/output/pext_future_source.txt examples/output/pext_future_source.pdf 50 10000
