# NLRexpress

Content summary:
* [General Info](#general-info)
* [Project status](#project-status)
* [Installation](#installation)
* [Usage](#usage)
* [Cite NLRexpress](#citing-nlrexpress)
* [References](#references)


## General Info
NLRexpress is an open-source bundle of 17 machine learning-based predictors trained for detecting sequence motifs specific to plant NLR-associated domains and is organised in 4 distinct modules:

* __CCexpress__: for an extended version of the EDVID motif: "_RDbbbDbEDVID_" (b - hydrophobic residue), which covers the 3rd alpha helice of the CC domain.
* __TIRexpress__: comprises 6 predictors for the following conserved motifs found in plant TIR domains : 
    * $\beta$ A: FLSFRGEDTR
    * $\alpha$ A: KNFTSHL
    * $\betaC$ : SRISIVVF
    * $\alphaC$: WCLDEL
    * $\betaD\alphaD1$: VLPVFYDVDPSDVRKQ
    * $\alphaD3$: WREALTEVANLSG

* __NBSexpress__: comprises 11 predictors for the following conserved motifs specific to plant NBS/NB-ARC domains : 
    * VG:            bbGRx
    * P-loop:        GbGGbGKTT
    * RNBS-A:        FDbrhWhshs
    * Walker-B:      KRbbbbDD
    * RNBS-B:        KbbbTTR
    * RNBS-C:        LseeeSWeLF
    * RNBS-D:        CFLYCSLFP
    * GLPL:          GLPLA
    * MHD:           bHD

* __LRRexpress__: for detecting "LxxLxL" motifs in Leucine-rich repeat (LRR) domains.

__Not limited to plant NLRs, the LRRexpress module was developed to be also used for screening other LRR-containing protein classes regardless of their taxonomy.__

The LRRexpress module is a much faster complementary alternative to LRRpredictor[1]. While LRRpredictor consists of a bundle of 8 voting estimators based on secondary structure predictions and sequence variability profiles generated on a large universal sequence database (Uniprot20), LRRexpress does not require the computation of secondary structure predictions and utilises a miniature targeted sequence database using JackHMMER [2][3] which significantly improves the speed of the pipeline with an average of less than 40s per NLR sequence (though average times might vary depending on computer specs).

For a limited number of protein sequences (under 100 seq), predictions can be run online on our webserver at https://nlrexpress.biochim.ro/ 



## Installation & initial setup 

NLRexpress was tested only on Ubuntu 18.04, but older/newer distributions should be compatible. 

#### 1. Cloning the project:
```
git clone https://github.com/eliza-m/NLRexpress
```

#### 2. Create conda environment
```
cd NLRexpress
conda env create -f environment.yml
conda activate nlrexpress
```

If conda (anaconda/miniconda) is not installed, installation steps can be found at https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html
	
#### 3. Install HMMER v3.3 (Nov 2019)

The HMMER version 3.3 can be installed in Ubuntu 18 via: 
```
sudo apt-get install hmmer
```

However, this might be superseded in the future by up-coming releases. The HMMER version archive can be found at : http://hmmer.org/download.html


#### 4. Download predictor models
```
  wget https://nlrexpress.biochim.ro/datasets/models.tar.gz
  tar -xf models.tar.gz
```
	
Now the installation is complete and NLRexpress ready for use. 

  
  
## Usage

The input file can contain multiple protein sequences in FASTA format. There is no constraint for the number of characters per line. The name of the query protein will be defined as the string beginning after the first '>' symbol until the first occurrence of a space character.

An example is provided with the repository. (Please note that the NLRexpress requires the conda environment to be active) :

While running a status log will be printed on the screen.

```
user@server:/home/user/git//NLRexpress$ conda activate nlrexpress

(nlrexpress) user@server:/home/user/git//NLRexpress$ python nlrexpress.py --input sample/input/zar1_rpp1.fa --outdir sample/output_ref/ --module all

22/08/2022 15:57:49:    ############ NLRexpress started ############
22/08/2022 15:57:49:    Input FASTA: sample/input/zar1_rpp1.fa
22/08/2022 15:57:49:    Checking FASTA file - started
22/08/2022 15:57:49:    Checking FASTA file - done
22/08/2022 15:57:49:    Running JackHMMER - started
22/08/2022 15:58:48:    Running JackHMMER - done
22/08/2022 15:58:48:    Preparing features: Parsing HMM profile - started
22/08/2022 15:58:48:    Running CCexpress : started
22/08/2022 15:58:49:    Preparing features: NN input for motif extEDVID started
22/08/2022 15:58:49:    Preparing features: NN input for motif extEDVID done
22/08/2022 15:58:49:    Running CCexpress : Running extEDVID predictor: ...done
22/08/2022 15:58:49:    Running CCexpress : done
22/08/2022 15:58:49:    Running TIRexpress : started
22/08/2022 15:58:49:    Preparing features: NN input for motif bA started
22/08/2022 15:58:49:    Preparing features: NN input for motif bA done
22/08/2022 15:58:49:    Running TIRexpress : Running bA predictor: ...done
22/08/2022 15:58:49:    Preparing features: NN input for motif aA started
22/08/2022 15:58:49:    Preparing features: NN input for motif aA done
22/08/2022 15:58:49:    Running TIRexpress : Running aA predictor: ...done
22/08/2022 15:58:49:    Preparing features: NN input for motif bC started
22/08/2022 15:58:50:    Preparing features: NN input for motif bC done
22/08/2022 15:58:50:    Running TIRexpress : Running bC predictor: ...done
22/08/2022 15:58:50:    Preparing features: NN input for motif aC started
22/08/2022 15:58:50:    Preparing features: NN input for motif aC done
22/08/2022 15:58:50:    Running TIRexpress : Running aC predictor: ...done
22/08/2022 15:58:50:    Preparing features: NN input for motif bDaD1 started
22/08/2022 15:58:50:    Preparing features: NN input for motif bDaD1 done
22/08/2022 15:58:50:    Running TIRexpress : Running bDaD1 predictor: ...done
22/08/2022 15:58:50:    Preparing features: NN input for motif aD3 started
22/08/2022 15:58:50:    Preparing features: NN input for motif aD3 done
22/08/2022 15:58:50:    Running TIRexpress : Running aD3 predictor: ...done
22/08/2022 15:58:50:    Running TIRexpress ...done
22/08/2022 15:58:50:    Running NBSexpress : started
22/08/2022 15:58:50:    Preparing features: NN input for motif VG started
22/08/2022 15:58:50:    Preparing features: NN input for motif VG done
22/08/2022 15:58:51:    Running NBSexpress : Running VG predictor: ...done
22/08/2022 15:58:51:    Preparing features: NN input for motif P-loop started
22/08/2022 15:58:51:    Preparing features: NN input for motif P-loop done
22/08/2022 15:58:52:    Running NBSexpress : Running P-loop predictor: ...done
22/08/2022 15:58:52:    Preparing features: NN input for motif RNSB-A started
22/08/2022 15:58:52:    Preparing features: NN input for motif RNSB-A done
22/08/2022 15:58:52:    Running NBSexpress : Running RNSB-A predictor: ...done
22/08/2022 15:58:52:    Preparing features: NN input for motif RNSB-B started
22/08/2022 15:58:52:    Preparing features: NN input for motif RNSB-B done
22/08/2022 15:58:53:    Running NBSexpress : Running RNSB-B predictor: ...done
22/08/2022 15:58:53:    Preparing features: NN input for motif RNSB-C started
22/08/2022 15:58:53:    Preparing features: NN input for motif RNSB-C done
22/08/2022 15:58:54:    Running NBSexpress : Running RNSB-C predictor: ...done
22/08/2022 15:58:54:    Preparing features: NN input for motif RNSB-D started
22/08/2022 15:58:54:    Preparing features: NN input for motif RNSB-D done
22/08/2022 15:58:54:    Running NBSexpress : Running RNSB-D predictor: ...done
22/08/2022 15:58:54:    Preparing features: NN input for motif Walker-B started
22/08/2022 15:58:54:    Preparing features: NN input for motif Walker-B done
22/08/2022 15:58:55:    Running NBSexpress : Running Walker-B predictor: ...done
22/08/2022 15:58:55:    Preparing features: NN input for motif GLPL started
22/08/2022 15:58:55:    Preparing features: NN input for motif GLPL done
22/08/2022 15:58:55:    Running NBSexpress : Running GLPL predictor: ...done
22/08/2022 15:58:55:    Preparing features: NN input for motif MHD started
22/08/2022 15:58:55:    Preparing features: NN input for motif MHD done
22/08/2022 15:58:56:    Running NBSexpress : Running MHD predictor: ...done
22/08/2022 15:58:56:    Running NBSexpress ...done
22/08/2022 15:58:56:    Running LRRexpress ...started
22/08/2022 15:58:56:    Preparing features: NN input for motif LxxLxL started
22/08/2022 15:58:56:    Preparing features: NN input for motif LxxLxL done
22/08/2022 15:58:56:    Running LRRexpress : Running LxxLxL predictor: ...done
22/08/2022 15:58:56:    Running LRRexpress ...done
22/08/2022 15:58:56:    Printing final results...started
22/08/2022 15:58:56:    Printing final results...done
22/08/2022 15:58:56:    ############ NLRexpress finished ############

```

Customized parameters are described in help mode:

```
python nlrexpress.py --help
Usage: nlrexpress.py [OPTIONS]

  Predict NLR-related motifs

Options:
  --input TEXT              Input FASTA file   [required]
  --outdir TEXT             Output folder  [required]
  --module TEXT             Predifined prediction modules :
                              - cc  :  CCexpress contains motif predictors for the CC domain:
                                                        extEDVID: rdhhhdhEDVID

                              - tir  :  TIRexpress contains motif predictors for the TIR domain:
                                                        bA: FLSFRGEDTR
                                                        aA: KNFTSHL
                                                        bC: SRISIVVF
                                                        aC: WCLDEL
                                                        bDaD1: VLPVFYDVDPSDVRKQ
                                                        aD3: WREALTEVANLSG


                              - nbs :  NBSexpress contains motif predictors for the NBS/NBARC domain:
                                                        VG:            bbGRx
                                                        P-loop:        GbGGbGKTT
                                                        RNBS-A:        FDbrhWhshs
                                                        Walker-B:      KRbbbbDD
                                                        RNBS-B:        KbbbTTR
                                                        RNBS-C:        LseeeSWeLF
                                                        RNBS-D:        CFLYCSLFP
                                                        GLPL:          GLPLA
                                                        MHD:           bHD

                              - lrr  :  LRRexpress contains motif predictors for the NBS/NBARC domain:
                                                        LRR pattern:   LxxLxL

                              - all  :  All modules above: CC, TIR, NBS and LRR express            [required]
                              
  --outformat TEXT          Output format layout :
                              - short          : Only the hits with more than 20% probability will be printed ( one line per motif )
                              - long           : All predicted residues will be printed ( one line per residue )
                              - all [default]  : Both short and long output formats are generated
                              
  --cpunum INTEGER          Set the number of CPU threads to be used.
  
  --writeinputfile BOOLEAN  Write input file to be used in case of job intreruption.
  
  --help                    Show this message and exit.
  
```




Prediction results will be found within the provided output directory path:

```
NLRexpress$ ls -l sample/output/
ls -l sample/output_ref/
-rw-rw-r-- 1 eliza eliza 965535 Aug 22 15:58 zar1_rpp1-1.hmm
-rw-rw-r-- 1 eliza eliza 965525 Aug 22 15:58 zar1_rpp1-2.hmm
-rw-rw-r-- 1 eliza eliza   2100 Aug 22 15:57 zar1_rpp1.fasta_proc
-rw-rw-r-- 1 eliza eliza   5360 Aug 22 15:58 zar1_rpp1.log
-rw-rw-r-- 1 eliza eliza 447768 Aug 22 15:58 zar1_rpp1.long.output.txt
-rw-rw-r-- 1 eliza eliza   7687 Aug 22 15:58 zar1_rpp1.short.output.txt
```

A sample output for `zar1_rpp1.fa` example can be found in `sample/output_ref/` for comparison purposes.

Two types of output formats are provided:  a summarising and a detailed version.

The summarizing (short) version `filename.short.output.txt` which contains from right to left: sequence name (from the FASTA header), the residue ID at which the predicted motif starts, the motif name, estimated probability, 5 positions upstream of the motif, the motif sequence and 5 positions downstream. Shown are only potential motifs yielding a probability estimate above 20%.

```
#ProtName                                           ResId           Motif    Proba   | -5pos |        MotifSeq | +5pos |
Q38834_ZAR1                                            66        extEDVID    99.74   | LVADL |    RELVYEAEDILV | DCQLA |
Q38834_ZAR1                                           160              VG    99.24   | YDHTQ |           VVGLE | GDKRK |
Q38834_ZAR1                                           189          P-loop   100.00   | IMAFV |       GMGGLGKTT | IAQEV |
Q38834_ZAR1                                           212          RNSB-A    99.81   | EIEHR |      FERRIWVSVS | QTFTE |
Q38834_ZAR1                                           260        Walker-B    97.30   | QYLLG |        KRYLIVMD | DVWDK |
Q38834_ZAR1                                           291          RNSB-B    98.46   | RGQGG |         SVIVTTR | SESVA |
Q38834_ZAR1                                           318          RNSB-C    99.94   | HRPEL |      LSPDNSWLLF | CNVAF |
Q38834_ZAR1                                           357            GLPL    99.98   | VTKCK |           GLPLT | IKAVG |
Q38834_ZAR1                                           418          RNSB-D    98.75   | SHLKS |       CILTLSLYP | EDCVI |
Q38834_ZAR1                                           487             MHD    99.65   | IITCK |             IHD | MVRDL |
Q38834_ZAR1                                           512          LxxLxL    93.98   | PEGLN |          CRHLGI | SGNFD |
Q38834_ZAR1                                           561          LxxLxL    99.73   | TDCKY |          LRVLDI | SKSIF |
Q38834_ZAR1                                           588          LxxLxL    99.93   | ASLQH |          LACLSL | SNTHP |
Q38834_ZAR1                                           612          LxxLxL    99.95   | EDLHN |          LQILDA | SYCQN |
Q38834_ZAR1                                           636          LxxLxL    99.90   | VLFKK |          LLVLDM | TNCGS |
Q38834_ZAR1                                           660          LxxLxL    54.88   | GSLVK |          LEVLLG | FKPAR |
Q38834_ZAR1                                           686          LxxLxL    99.87   | KNLTN |          LRKLGL | SLTRG |
...
```

The detailed long format `filename.long.output.txt`, contains one row per sequence residue and the corresponding probability estimate for each of the motif predictors. The first and last positions of the sequences do not have an estimate due to the fact that the predictors require 5 positions up- and down-stream the predicted position.

```
#ProtName                                           resid  aa   extEDVID              bA              aA              bC              a   bDaD1      aD3              VG          P-loop          RNSB-A          RNSB-B          RNSB-C          RNSB-D        Walker-B           GLPL      MHD          LxxLxL
Q38834_ZAR1                                             1   M          -               -               -               -
      -        -               -               -               -               -               -               -               -              -        -               -
Q38834_ZAR1                                             2   V          -               -               -               -
      -        -               -               -               -               -               -               -               -              -        -               -
Q38834_ZAR1                                             3   D          -               -               -               -
      -        -               -               -               -               -               -               -               -              -        -               -
Q38834_ZAR1                                             4   A          -               -               -               -
      -        -               -               -               -               -               -               -               -              -        -               -
Q38834_ZAR1                                             5   V          -               -               -               -
      -        -               -               -               -               -               -               -               -              -        -               -
Q38834_ZAR1                                             6   V       0.00            0.00            0.00            0.00            0.0    0.00     0.00            0.04            0.00            0.00            0.02            0.00            0.00            0.00           0.00     0.09            0.00
Q38834_ZAR1                                             7   T       0.00            0.00            0.01            0.00            0.0    0.00     0.00            0.02            0.00            0.00            0.01            0.00            0.00            0.00           0.00     0.00            0.00
Q38834_ZAR1                                             8   V       0.00            0.00            0.00            0.00            0.0    0.00     0.00            0.02            0.00            0.00            0.00            0.00            0.00            0.00           0.00     0.00            0.00
Q38834_ZAR1                                             9   F       0.00            0.00            0.00            0.00            0.0    0.00     0.00            0.03            0.00            0.00            0.00            0.00            0.00            0.00           0.00     0.00            0.00
Q38834_ZAR1                                            10   L       0.00            0.00            0.00            0.00            0.0    0.00     0.00            0.01            0.00            0.00            0.00            0.00            0.00            0.00           0.00     0.07            0.00
Q38834_ZAR1                                            11   E       0.00            0.00            0.00            0.00            0.0    0.00     0.00            0.00            0.00            0.00            0.00            0.00            0.00            0.00           0.00     0.00            0.00
Q38834_ZAR1                                            12   K       0.00            0.00            0.00            0.00            0.0    0.00     0.00            0.00            0.00            0.00            0.02            0.00            0.00            0.00           0.00     0.00            0.00
Q38834_ZAR1                                            13   T       0.00            0.00            0.00            0.00            0.0    0.00     0.00            0.02            0.00            0.00            0.01            0.00            0.00            0.00           0.00     0.02            0.10
Q38834_ZAR1                                            14   L       0.10            0.00            0.00            0.00            0.0    0.00     0.00            0.04            0.00            0.00            0.02            0.00            0.00            0.00           0.00     0.00            0.08
Q38834_ZAR1                                            15   N       0.00            0.00            0.00            0.00            0.0    0.00     0.00            0.00            0.00            0.00            0.00            0.00            0.00            0.00           0.00     0.00            0.00
Q38834_ZAR1                                            16   I       0.00            0.00            0.00            0.00            0.0    0.00     0.00            0.18            0.00            0.00            0.00            0.00            0.00            0.00           0.00     0.01            0.01
Q38834_ZAR1                                            17   L       0.00            0.00            0.00            0.00            0.0    0.00     0.00            0.00            0.00            0.00            0.00            0.00            0.00            0.00           0.00     0.00            0.00
Q38834_ZAR1                                            18   E       0.00            0.00            0.00            0.00            0.0    0.00     0.00            0.00            0.00            0.00            0.00            0.00            0.00            0.00           0.00     0.00            0.00
Q38834_ZAR1                                            19   E       0.00            0.00            0.00            0.00            0.0    0.00     0.00            0.00            0.00            0.00            0.00            0.00            0.00            0.00           0.00     0.00            0.00
```



## Citing LRRpredictor

If you use LRRpredictor please cite:

[Eliza C. Martin, Laurentiu Spiridon, Aska Goverse, Andrei-Jose Petrescu. NLRexpress - a bundle of machine learning motif predictors - reveals motif stability underlying plant NLRs diversity. Frontiers in Plant Science 2022.](https://www.frontiersin.org/articles/10.3389/fpls.2022.975888/abstract)


## References

[1] Martin EC, Sukarta OCA, Spiridon L, Grigore LG, Constantinescu V, Tacutu R, Goverse A, Petrescu AJ. LRRpredictor-A New LRR Motif Detection Method for Irregular Motifs of Plant NLR Proteins Using an Ensemble of Classifiers. Genes (Basel). 2020 Mar 8;11(3):286. doi: 10.3390/genes11030286.

[2] Jaina Mistry, Robert D. Finn, Sean R. Eddy, Alex Bateman, Marco Punta, Challenges in homology search: HMMER3 and convergent evolution of coiled-coil regions, Nucleic Acids Research, Volume 41, Issue 12, 1 July 2013, Page e121, https://doi.org/10.1093/nar/gkt263.

[3] Robert D. Finn, Jody Clements, Sean R. Eddy, HMMER web server: interactive sequence similarity searching, Nucleic Acids Research, Volume 39, Issue suppl_2, 1 July 2011, Pages W29–W37, https://doi.org/10.1093/nar/gkr367

