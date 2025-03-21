# Codon_optimizer by Sebastian Einhauser
A simple yet powerful console based tool for codon optimization of DNA sequences.
## Functionality
The tool reads in DNA or AA sequences provided as text.
After detection of possible reading frames, and user input, to define the wantend reading frame, multiple Codon Usage profiles are available. On top, multiple (de)optimization algorithms are available.
After base optimization, several other optimization options like splice site avoidance, slippery site avoidance, custom site avoidance or one-to stop codon avoidance are available.

## Codon Usage Profiles
There are 5 built-in Codon Usage profiles available:
* Human
* E.Coli
* Mouse
* Yeast(Saccharomyces)
* Insect
  
## Optimization algorithms
There are seven algoritms are available:
#### Optimization:
* Most Frequent: using only the most prevalent codon for that aminoacid.
* Probability Distribution: Each codon is sampled from the according probability distribution of all available codons for that aa.
* Enforced Distribution: Tries too mimic the chosen organism codon profile as close as possible.

#### Deoptimization:
* Least Frequent: using only the least prevalent codon for that aminoacid.
* Inverted Probability Distribution: Each codon is sampled from the inverse probability distribution of all available codons for that aa making low frequency high and vice versa.
* Enforced inverted Distribution: Tries too mimic the chosen organism inverse codon profile as close as possible.
#### Special:
* GeneOptimizer: Inspired By Raab et al. 2009, Applying codon usage, GC optimization and repetition avoidance in a sliding window approach. Calculating all possible sequences in this window and score them with a scoring function. User can choose target GC content.
  
## Additional functionality
There are 5 additional built in Functions:
* Codon Duplet avoidance: Avoids using two times the same codon in a row.
* Slippery Site Avoidance: Avoids 4x the same base in a row to prevent frame shifts.
* Splice Site Avoidance:
- Cryptic Sites: Avoids a 24base window of known cryptic splice sites (CPU intensive, 1 Million windows will be scanned)
- Consensus Sites: Avoids the consensus eukariotic splice site.
* One-To-Stop Codons:
  - avoid: avoids ots codons, so codons only taking a single point mutation to become a stop.
  - implement: Mutates all codons to a OTS when possible.
* Custom Sites: Allows the user to input custom sequence motifs to avoid, e.g. restriciton sites.
## Additional Notes:
Additional functions should be ran at least 2 times as mutations introduced by the first run might result in the introduction of novel unwantend sites. This is true across functions but also within one function.
E.g. Splice-Site-avoidance may introduce slippery sites.
E.g. Slippery Site Avoidance may erase one Slippery site but accidently introduce a new one a few bases before.
