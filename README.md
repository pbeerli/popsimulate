# popsimulate
Simulator of population genetics data allowing for multiple populations and immigration

To generate simulated datasets do this:

  migtree < parmfile.stepstone | migdata 

It will print some auxillary stuff to the screen but will 
produce a file called infile with the sequences in it.

migdata is not interactive and gets all its stuff through the migtree
output what that is you can check by doing

migtree < parmfile.stepstone > gugus

and then look at gugus

Here an example parmfile for the simlator
I explain each line after #
You will need to remove the #comments to make the file work
        
  S Nogamma          # {S: datatype [S, A, M = seq, allelic, microsat]} 
                     # {Nogamma:  mutation rate is NOT gamma deviated}
                     # { if you want to use Gamma let me know [->more explan}
  4 3 9839017        # number_of_populations number_of_loci random_number_seed
                     # the random number seed should be of the type 4n+1
   10 10 10 10       # sample size for each population
  100 100 100 100    # population size for each population, remember this
                     # this gets multiplied by and 4 -> theta=4 N mu
  - 0 0 0            # immigration matrix - on diagonal 
  0.00125 - 0 0      # this example has no migration except from
  0 0.00125 - 0      #  1->2->3->4
  0 0 0.00125 -      #
  0.0000125          # mutation rate per site
  1000 2.0           # number_of_sites t/t_ratio
