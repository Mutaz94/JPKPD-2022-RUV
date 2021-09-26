Sat Sep 25 01:06:56 CDT 2021
$PROB template control stream
;-----------------------------------------------------------------------
; Project: 	Investigating the contribution of residual unexplained
; 	   	variability in nonlinear mixed-effect approach
; Model: 	Two-compartment model with linear elimination
; Estim:	First-order conditional est. with interaction
; Author: 	Mutaz M. Jaber <jaber038@umn.edu>
; Date created: 9/7/2021
; Date modified: 9/7/2021
;-----------------------------------------------------------------------
$INPUT ID TIME DV AMT MDV EVID
$DATA ../../../../data/int/SL2/dat29.csv ignore=@
$SUBR ADVAN4 TRANS4
$EST MET=1 NOABORT MAX=10000 PRINT=5 INTER
$PK
ET1 = EXP(ETA(1)*THETA(6))
ET2 = EXP(ETA(2)*THETA(7))
ET3 = EXP(ETA(3)*THETA(8))
ET4 = EXP(ETA(4)*THETA(9))
ET5 = EXP(ETA(5)*THETA(10))

CL = 5.0 * THETA(1) * ET1
V2 = 35  * THETA(2) * ET2
Q  = 50  * THETA(3) * ET3
V3 = 50  * THETA(4) * ET4
KA = 0.7 * THETA(5) * ET5
SC = V2
$ERROR
CVERR = 0.05
W = THETA(11)*F*CVERR

Y 	= F + W*ERR(1)

$THETA
(0,1) ; CL
(0,1) ; V2
(0,1) ; Q
(0,1) ; V3
(0,1) ; KA
(0,1) ; IIVCL
(0,1) ; IIVV2
(0,1) ; IIVQ
(0,1) ; IIVV3
(0,1) ; IIVKA
(0,1) ; CVPropErr

$OMEGA  0.09 FIX 0.09 FIX 0.09 FIX 0.09 FIX 0.09 FIX
$SIGMA  1 FIX ;        [P]
$COVARIANCE UNCOND

NM-TRAN MESSAGES
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.

License Registered to: University of Minnesota
Expiration Date:    14 APR 2022
Current Date:       25 SEP 2021
Days until program expires : 204
1NONLINEAR MIXED EFFECTS MODEL PROGRAM (NONMEM) VERSION 7.5.0
 ORIGINALLY DEVELOPED BY STUART BEAL, LEWIS SHEINER, AND ALISON BOECKMANN
 CURRENT DEVELOPERS ARE ROBERT BAUER, ICON DEVELOPMENT SOLUTIONS,
 AND ALISON BOECKMANN. IMPLEMENTATION, EFFICIENCY, AND STANDARDIZATION
 PERFORMED BY NOUS INFOSYSTEMS.

 PROBLEM NO.:         1
 template control stream
0DATA CHECKOUT RUN:              NO
 DATA SET LOCATED ON UNIT NO.:    2
 THIS UNIT TO BE REWOUND:        NO
 NO. OF DATA RECS IN DATA SET:      997
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (2E4.0,E20.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      897
 TOT. NO. OF INDIVIDUALS:      100
0LENGTH OF THETA:  11
0DEFAULT THETA BOUNDARY TEST OMITTED:    NO
0OMEGA HAS SIMPLE DIAGONAL FORM WITH DIMENSION:   5
0DEFAULT OMEGA BOUNDARY TEST OMITTED:    NO
0SIGMA HAS SIMPLE DIAGONAL FORM WITH DIMENSION:   1
0DEFAULT SIGMA BOUNDARY TEST OMITTED:    NO
0INITIAL ESTIMATE OF THETA:
 LOWER BOUND    INITIAL EST    UPPER BOUND
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
0INITIAL ESTIMATE OF OMEGA:
 0.9000E-01
 0.0000E+00   0.9000E-01
 0.0000E+00   0.0000E+00   0.9000E-01
 0.0000E+00   0.0000E+00   0.0000E+00   0.9000E-01
 0.0000E+00   0.0000E+00   0.0000E+00   0.0000E+00   0.9000E-01
0OMEGA CONSTRAINED TO BE THIS INITIAL ESTIMATE
0INITIAL ESTIMATE OF SIGMA:
 0.1000E+01
0SIGMA CONSTRAINED TO BE THIS INITIAL ESTIMATE
0COVARIANCE STEP OMITTED:        NO
 EIGENVLS. PRINTED:              NO
 SPECIAL COMPUTATION:            NO
 COMPRESSED FORMAT:              NO
 GRADIENT METHOD USED:     NOSLOW
 SIGDIGITS ETAHAT (SIGLO):                  -1
 SIGDIGITS GRADIENTS (SIGL):                -1
 EXCLUDE COV FOR FOCE (NOFCOV):              NO
 Cholesky Transposition of R Matrix (CHOLROFF):0
 KNUTHSUMOFF:                                -1
 RESUME COV ANALYSIS (RESUME):               NO
 SIR SAMPLE SIZE (SIRSAMPLE):
 NON-LINEARLY TRANSFORM THETAS DURING COV (THBND): 1
 PRECONDTIONING CYCLES (PRECOND):        0
 PRECONDTIONING TYPES (PRECONDS):        TOS
 FORCED PRECONDTIONING CYCLES (PFCOND):0
 PRECONDTIONING TYPE (PRETYPE):        0
 FORCED POS. DEFINITE SETTING DURING PRECONDITIONING: (FPOSDEF):0
 SIMPLE POS. DEFINITE SETTING: (POSDEF):-1
1DOUBLE PRECISION PREDPP VERSION 7.5.0

 TWO COMPARTMENT MODEL WITH FIRST-ORDER ABSORPTION (ADVAN4)
0MAXIMUM NO. OF BASIC PK PARAMETERS:   5
0BASIC PK PARAMETERS (AFTER TRANSLATION):
   BASIC PK PARAMETER NO.  1: ELIMINATION RATE (K)
   BASIC PK PARAMETER NO.  2: CENTRAL-TO-PERIPH. RATE (K23)
   BASIC PK PARAMETER NO.  3: PERIPH.-TO-CENTRAL RATE (K32)
   BASIC PK PARAMETER NO.  5: ABSORPTION RATE (KA)
 TRANSLATOR WILL CONVERT PARAMETERS
 CL, V2, Q, V3 TO K, K23, K32 (TRANS4)
0COMPARTMENT ATTRIBUTES
 COMPT. NO.   FUNCTION   INITIAL    ON/OFF      DOSE      DEFAULT    DEFAULT
                         STATUS     ALLOWED    ALLOWED    FOR DOSE   FOR OBS.
    1         DEPOT        OFF        YES        YES        YES        NO
    2         CENTRAL      ON         NO         YES        NO         YES
    3         PERIPH.      ON         NO         YES        NO         NO
    4         OUTPUT       OFF        YES        NO         NO         NO
1
 ADDITIONAL PK PARAMETERS - ASSIGNMENT OF ROWS IN GG
 COMPT. NO.                             INDICES
              SCALE      BIOAVAIL.   ZERO-ORDER  ZERO-ORDER  ABSORB
                         FRACTION    RATE        DURATION    LAG
    1            *           *           *           *           *
    2            6           *           *           *           *
    3            *           *           *           *           *
    4            *           -           -           -           -
             - PARAMETER IS NOT ALLOWED FOR THIS MODEL
             * PARAMETER IS NOT SUPPLIED BY PK SUBROUTINE;
               WILL DEFAULT TO ONE IF APPLICABLE
0DATA ITEM INDICES USED BY PRED ARE:
   EVENT ID DATA ITEM IS DATA ITEM NO.:      6
   TIME DATA ITEM IS DATA ITEM NO.:          2
   DOSE AMOUNT DATA ITEM IS DATA ITEM NO.:   4

0PK SUBROUTINE CALLED WITH EVERY EVENT RECORD.
 PK SUBROUTINE NOT CALLED AT NONEVENT (ADDITIONAL OR LAGGED) DOSE TIMES.
0ERROR SUBROUTINE CALLED WITH EVERY EVENT RECORD.
1


 #TBLN:      1
 #METH: First Order Conditional Estimation with Interaction

 ESTIMATION STEP OMITTED:                 NO
 ANALYSIS TYPE:                           POPULATION
 NUMBER OF SADDLE POINT RESET ITERATIONS:      0
 GRADIENT METHOD USED:               NOSLOW
 CONDITIONAL ESTIMATES USED:              YES
 CENTERED ETA:                            NO
 EPS-ETA INTERACTION:                     YES
 LAPLACIAN OBJ. FUNC.:                    NO
 NO. OF FUNCT. EVALS. ALLOWED:            10000
 NO. OF SIG. FIGURES REQUIRED:            3
 INTERMEDIATE PRINTOUT:                   YES
 ESTIMATE OUTPUT TO MSF:                  NO
 ABORT WITH PRED EXIT CODE 1:             NO
 IND. OBJ. FUNC. VALUES SORTED:           NO
 NUMERICAL DERIVATIVE
       FILE REQUEST (NUMDER):               NONE
 MAP (ETAHAT) ESTIMATION METHOD (OPTMAP):   0
 ETA HESSIAN EVALUATION METHOD (ETADER):    0
 INITIAL ETA FOR MAP ESTIMATION (MCETA):    0
 SIGDIGITS FOR MAP ESTIMATION (SIGLO):      100
 GRADIENT SIGDIGITS OF
       FIXED EFFECTS PARAMETERS (SIGL):     100
 NOPRIOR SETTING (NOPRIOR):                 0
 NOCOV SETTING (NOCOV):                     OFF
 DERCONT SETTING (DERCONT):                 OFF
 FINAL ETA RE-EVALUATION (FNLETA):          1
 EXCLUDE NON-INFLUENTIAL (NON-INFL.) ETAS
       IN SHRINKAGE (ETASTYPE):             NO
 NON-INFL. ETA CORRECTION (NONINFETA):      0
 RAW OUTPUT FILE (FILE): m29.ext
 EXCLUDE TITLE (NOTITLE):                   NO
 EXCLUDE COLUMN LABELS (NOLABEL):           NO
 FORMAT FOR ADDITIONAL FILES (FORMAT):      S1PE12.5
 PARAMETER ORDER FOR OUTPUTS (ORDER):       TSOL
 KNUTHSUMOFF:                               0
 INCLUDE LNTWOPI:                           NO
 INCLUDE CONSTANT TERM TO PRIOR (PRIORC):   NO
 INCLUDE CONSTANT TERM TO OMEGA (ETA) (OLNTWOPI):NO
 ADDITIONAL CONVERGENCE TEST (CTYPE=4)?:    NO
 EM OR BAYESIAN METHOD USED:                 NONE


 THE FOLLOWING LABELS ARE EQUIVALENT
 PRED=PREDI
 RES=RESI
 WRES=WRESI
 IWRS=IWRESI
 IPRD=IPREDI
 IRS=IRESI

 MONITORING OF SEARCH:


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1459.37475787304        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   6.8262E+01 -7.7774E+01  6.6362E+01  5.5656E+01  1.3717E+02 -2.9871E+01 -1.4970E+02 -1.7109E+02 -1.2152E+02 -4.2494E+01
            -4.3685E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2855.30382624223        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.0218E+00  1.3326E+00  1.0649E+00  7.6710E-01  1.2014E+00  1.1180E+00  1.4653E+00  9.5933E-01  1.1683E+00  1.1690E+00
             1.9966E+00
 PARAMETER:  1.2159E-01  3.8712E-01  1.6284E-01 -1.6514E-01  2.8348E-01  2.1158E-01  4.8208E-01  5.8476E-02  2.5555E-01  2.5617E-01
             7.9145E-01
 GRADIENT:   5.6184E+01 -3.7047E+01 -1.2692E+01 -5.0379E+01 -1.7132E+00  1.5993E+01  3.2373E+01 -4.6385E+00 -6.1875E+00 -8.7917E+00
            -3.7282E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2860.24065452675        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  9.5941E-01  1.7688E+00  3.4877E+00  6.6286E-01  2.0464E+00  1.2423E+00  1.3497E+00  4.0134E-01  1.4307E+00  1.6999E+00
             2.0956E+00
 PARAMETER:  5.8562E-02  6.7031E-01  1.3492E+00 -3.1119E-01  8.1611E-01  3.1696E-01  3.9986E-01 -8.1294E-01  4.5818E-01  6.3057E-01
             8.3986E-01
 GRADIENT:  -4.6447E+01  6.8050E+01 -9.7780E+00  4.5627E+01  9.0883E+01  4.8408E+01  4.6162E+01 -2.5361E-01  1.6803E+01 -2.7273E+01
            -2.4950E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2893.36429704332        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      226
 NPARAMETR:  1.0023E+00  1.7233E+00  1.0663E+00  6.0151E-01  1.4475E+00  1.0453E+00  1.0320E+00  6.3083E-02  1.2835E+00  1.4812E+00
             2.3174E+00
 PARAMETER:  1.0227E-01  6.4426E-01  1.6421E-01 -4.0831E-01  4.6985E-01  1.4429E-01  1.3150E-01 -2.6633E+00  3.4958E-01  4.9287E-01
             9.4044E-01
 GRADIENT:   8.2198E+00  7.5157E+00  2.6994E+00  7.6058E-01 -1.1861E+01 -7.8100E+00 -1.6420E+00 -1.7274E-02  3.2124E-01  4.0581E-01
             5.3621E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2894.22368302098        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      297
 NPARAMETR:  9.9789E-01  1.7758E+00  1.0477E+00  5.6590E-01  1.5022E+00  1.0599E+00  1.0146E+00  5.8313E-02  1.3205E+00  1.5136E+00
             2.3143E+00
 PARAMETER:  9.7891E-02  6.7426E-01  1.4657E-01 -4.6934E-01  5.0695E-01  1.5820E-01  1.1452E-01 -2.7419E+00  3.7804E-01  5.1449E-01
             9.3909E-01
 GRADIENT:  -2.1329E-01  4.8569E+00  1.5372E+00  1.3997E+00 -3.2927E+00 -2.5323E+00 -1.3839E+00 -1.3592E-02  8.8039E-01  6.6173E-01
             4.6068E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2897.01424309485        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      478
 NPARAMETR:  1.0069E+00  2.1855E+00  6.2164E-01  3.2623E-01  1.7448E+00  1.0795E+00  9.0024E-01  3.5689E-02  1.7372E+00  1.6813E+00
             2.2989E+00
 PARAMETER:  1.0688E-01  8.8183E-01 -3.7540E-01 -1.0201E+00  6.5666E-01  1.7651E-01 -5.0914E-03 -3.2329E+00  6.5230E-01  6.1956E-01
             9.3244E-01
 GRADIENT:   7.7661E+00  2.3607E+01 -8.8217E-01  8.8709E+00  3.0299E+00  1.7651E+00 -8.8237E-01 -2.2095E-03 -6.9240E-01  1.4171E+00
            -4.1976E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2898.39906943325        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      654
 NPARAMETR:  1.0058E+00  2.4690E+00  2.1827E-01  1.2642E-01  1.9196E+00  1.0733E+00  8.2909E-01  1.8809E-02  2.9655E+00  1.7969E+00
             2.2837E+00
 PARAMETER:  1.0582E-01  1.0038E+00 -1.4220E+00 -1.9682E+00  7.5214E-01  1.7070E-01 -8.7426E-02 -3.8734E+00  1.1870E+00  6.8608E-01
             9.2578E-01
 GRADIENT:   6.1431E+00  1.5052E+01  1.1626E-01  5.5165E-01  2.0979E+00 -8.7568E-01  7.7854E-01  6.1252E-05 -1.2496E+00  1.1656E+00
            -3.6788E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2898.48681925142        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      829
 NPARAMETR:  1.0024E+00  2.4568E+00  2.1641E-01  1.2577E-01  1.9101E+00  1.0757E+00  8.2459E-01  1.7806E-02  3.0621E+00  1.7851E+00
             2.2858E+00
 PARAMETER:  1.0244E-01  9.9884E-01 -1.4306E+00 -1.9733E+00  7.4713E-01  1.7297E-01 -9.2870E-02 -3.9282E+00  1.2191E+00  6.7946E-01
             9.2670E-01
 GRADIENT:  -1.3787E-01  8.4461E-02 -4.9834E-02  7.8118E-02  1.3293E-01  1.4911E-02 -1.0386E-01  9.4284E-05 -4.6357E-04 -9.1147E-03
            -8.6522E-03

0ITERATION NO.:   39    OBJECTIVE VALUE:  -2898.48692743928        NO. OF FUNC. EVALS.: 129
 CUMULATIVE NO. OF FUNC. EVALS.:      958
 NPARAMETR:  1.0025E+00  2.4569E+00  2.1653E-01  1.2549E-01  1.9100E+00  1.0757E+00  8.2483E-01  1.7822E-02  3.0658E+00  1.7853E+00
             2.2858E+00
 PARAMETER:  1.0251E-01  9.9892E-01 -1.4300E+00 -1.9756E+00  7.4710E-01  1.7293E-01 -9.2583E-02 -3.9273E+00  1.2203E+00  6.7957E-01
             9.2670E-01
 GRADIENT:   3.1884E-03 -9.6049E-02 -5.5259E-04 -8.3268E-03 -2.0413E-02  3.5240E-03  5.2898E-03  9.4441E-05 -5.3723E-04  5.5180E-03
             1.3797E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      958
 NO. OF SIG. DIGITS IN FINAL EST.:  3.7

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.1399E-03 -1.8409E-02 -4.2292E-05  3.2064E-02 -1.8764E-02
 SE:             2.9518E-02  2.7719E-02  3.9430E-05  1.5382E-02  2.6744E-02
 N:                     100         100         100         100         100

 P VAL.:         9.6920E-01  5.0660E-01  2.8346E-01  3.7120E-02  4.8293E-01

 ETASHRINKSD(%)  1.1125E+00  7.1371E+00  9.9868E+01  4.8467E+01  1.0404E+01
 ETASHRINKVR(%)  2.2127E+00  1.3765E+01  1.0000E+02  7.3444E+01  1.9725E+01
 EBVSHRINKSD(%)  1.1665E+00  6.6446E+00  9.9858E+01  5.9286E+01  6.7492E+00
 EBVSHRINKVR(%)  2.3194E+00  1.2848E+01  1.0000E+02  8.3424E+01  1.3043E+01
 RELATIVEINF(%)  9.7661E+01  3.2705E+01  1.1628E-04  4.9258E+00  6.1564E+01
 EPSSHRINKSD(%)  1.7203E+01
 EPSSHRINKVR(%)  3.1446E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          897
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1648.5757285691827     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2898.4869274392781     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1249.9111988700954     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    22.47
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    13.73
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2898.487       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.00E+00  2.46E+00  2.17E-01  1.25E-01  1.91E+00  1.08E+00  8.25E-01  1.78E-02  3.07E+00  1.79E+00  2.29E+00
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


         ETA1      ETA2      ETA3      ETA4      ETA5     
 
 ETA1
+        9.00E-02
 
 ETA2
+        0.00E+00  9.00E-02
 
 ETA3
+        0.00E+00  0.00E+00  9.00E-02
 
 ETA4
+        0.00E+00  0.00E+00  0.00E+00  9.00E-02
 
 ETA5
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  9.00E-02
 


 SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS  ****


         EPS1     
 
 EPS1
+        1.00E+00
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


         ETA1      ETA2      ETA3      ETA4      ETA5     
 
 ETA1
+        3.00E-01
 
 ETA2
+        0.00E+00  3.00E-01
 
 ETA3
+        0.00E+00  0.00E+00  3.00E-01
 
 ETA4
+        0.00E+00  0.00E+00  0.00E+00  3.00E-01
 
 ETA5
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  3.00E-01
 


 SIGMA - CORR MATRIX FOR RANDOM EFFECTS - EPSILONS  ***


         EPS1     
 
 EPS1
+        1.00E+00
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                                     R MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        9.38E+02
 
 TH 2
+       -5.60E+00  2.40E+02
 
 TH 3
+        2.76E+00  2.55E+01  2.00E+02
 
 TH 4
+       -1.25E+01  2.76E+02 -3.50E+02  1.79E+03
 
 TH 5
+       -1.56E+00 -1.27E+01 -1.87E+01  1.02E+02  8.31E+01
 
 TH 6
+        3.72E+00 -2.43E+00  8.88E-01 -4.84E+00 -8.98E-01  1.63E+02
 
 TH 7
+       -4.35E+00  8.26E-02 -1.41E+01 -1.06E+01 -9.90E-01  8.28E-01  2.29E+02
 
 TH 8
+        5.66E+00  1.46E-01 -4.46E+00 -3.20E+00  1.20E-01 -2.44E+00  2.32E+00  1.96E+01
 
 TH 9
+       -2.97E-02 -1.35E+00 -3.62E+00  4.36E+01  3.74E-01 -1.26E-01  2.72E+00 -1.56E-02  3.26E+00
 
 TH10
+        2.32E-01 -3.11E+00  4.29E+04  3.03E+01 -6.21E+00  4.56E-01  8.72E-01  7.79E-01  1.04E+00  4.42E+01
 
 TH11
+       -1.25E+01 -1.02E+01 -4.34E-01 -3.96E+00  1.71E-01  3.22E+00  5.46E+00 -1.82E-01  1.58E+00  2.46E+00  2.23E+02
 
 OM11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM12
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         .........
 
 OM13
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... .........
 
 OM14
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... .........
 
 OM15
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... .........
 
 OM22
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... .........
 
 OM23
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... .........
 
 OM24
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... .........
 
1

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 OM25
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... .........
 
 OM33
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM34
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM35
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM44
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM45
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        .........
 
 OM55
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        ......... .........
 
 SG11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        ......... ......... .........
 
 Elapsed finaloutput time in seconds:     0.01
 #CPUT: Total CPU Time in Seconds,       36.312
Stop Time:
Sat Sep 25 01:07:34 CDT 2021
