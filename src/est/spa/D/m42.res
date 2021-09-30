Wed Sep 29 20:02:04 CDT 2021
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
$DATA ../../../../data/spa/D/dat42.csv ignore=@
$SUBR ADVAN4 TRANS4
$EST MET=1 NOABORT MAX=10000 PRINT=5 INTER NSIG=2
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

$OMEGA  (0.09 FIX)x5
$SIGMA  1 FIX ;        [P]
$COVARIANCE UNCOND

NM-TRAN MESSAGES
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.

License Registered to: University of Minnesota
Expiration Date:    14 APR 2022
Current Date:       29 SEP 2021
Days until program expires : 200
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
 NO. OF DATA RECS IN DATA SET:      500
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (E4.0,E3.0,E21.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      400
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
 NO. OF SIG. FIGURES REQUIRED:            2
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
 RAW OUTPUT FILE (FILE): m42.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   13745.2992035195        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.9288E+02  2.8727E+02 -3.1846E+01  3.8665E+02  1.2657E+02 -1.2097E+03 -7.1047E+02 -3.4419E+01 -8.8540E+02 -3.1143E+02
            -2.7313E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -662.132156913489        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.3609E+00  1.0469E+00  9.9771E-01  1.2479E+00  1.1541E+00  1.5833E+00  1.1665E+00  9.8439E-01  1.0210E+00  9.9357E-01
             1.5073E+01
 PARAMETER:  4.0815E-01  1.4583E-01  9.7706E-02  3.2149E-01  2.4330E-01  5.5950E-01  2.5399E-01  8.4270E-02  1.2081E-01  9.3554E-02
             2.8129E+00
 GRADIENT:   2.4406E+01 -1.7473E+01 -1.8752E+00 -3.5840E+01 -3.1130E-01  5.6816E+01  6.0727E+00  2.6832E+00  1.5637E+01  4.0903E+00
             2.0477E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -688.629160317944        NO. OF FUNC. EVALS.:  77
 CUMULATIVE NO. OF FUNC. EVALS.:      160
 NPARAMETR:  1.2253E+00  7.4935E-01  1.1332E+00  1.5093E+00  4.1093E+00  1.2672E+00  2.6117E+00  3.7850E-01  7.7536E-01  4.0114E+00
             1.3333E+01
 PARAMETER:  3.0322E-01 -1.8855E-01  2.2502E-01  5.1162E-01  1.5132E+00  3.3679E-01  1.0600E+00 -8.7154E-01 -1.5442E-01  1.4891E+00
             2.6902E+00
 GRADIENT:  -3.0808E+01  2.6406E+01  2.5550E+00  5.3955E+01 -4.7818E+00 -1.3605E+01  1.0351E+01  1.0401E-01  1.3468E+01 -6.0627E-01
             1.2292E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -715.134809870928        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      230
 NPARAMETR:  1.0445E+00  3.5176E-01  9.3085E-01  1.4307E+00  5.3309E+00  1.2261E+00  2.3380E+00  6.7635E-02  4.2016E-01  5.9926E+00
             1.0610E+01
 PARAMETER:  1.4358E-01 -9.4482E-01  2.8339E-02  4.5813E-01  1.7735E+00  3.0385E-01  9.4929E-01 -2.5936E+00 -7.6712E-01  1.8905E+00
             2.4618E+00
 GRADIENT:  -6.6479E+01  1.6262E+01  1.0916E+01  7.3802E+01 -9.2055E-01 -1.6946E+01  2.3647E+00  4.1120E-03  4.0107E+00 -3.4287E+00
            -1.5461E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -736.960260018804        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      302
 NPARAMETR:  9.1585E-01  1.5110E-01  1.7683E-01  1.0142E+00  1.0726E+01  1.2187E+00  2.9687E-01  1.0000E-02  3.6863E-02  2.2460E+00
             1.0274E+01
 PARAMETER:  1.2097E-02 -1.7898E+00 -1.6325E+00  1.1410E-01  2.4727E+00  2.9778E-01 -1.1145E+00 -1.2335E+01 -3.2006E+00  9.0914E-01
             2.4296E+00
 GRADIENT:   4.7988E+01  1.4130E+01 -5.8731E+01  1.4019E+02 -3.2373E+00 -4.1726E+01  5.6519E-01  0.0000E+00  3.5740E-02  5.7197E+00
            -8.5115E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -756.234719304655        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      375
 NPARAMETR:  7.0113E-01  6.1486E-02  9.5858E-02  6.8972E-01  1.3649E+01  1.3095E+00  3.1713E-02  1.0000E-02  1.0000E-02  1.5167E+00
             1.0687E+01
 PARAMETER: -2.5506E-01 -2.6890E+00 -2.2449E+00 -2.7146E-01  2.7137E+00  3.6962E-01 -3.3510E+00 -1.9209E+01 -6.3662E+00  5.1656E-01
             2.4690E+00
 GRADIENT:  -2.9659E+01 -7.2727E-01 -1.2172E+01  8.2993E+01  3.5945E-01  2.2032E+00  1.1528E-03  0.0000E+00  0.0000E+00  1.6948E-02
             8.7978E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -757.486024614176        NO. OF FUNC. EVALS.: 144
 CUMULATIVE NO. OF FUNC. EVALS.:      519
 NPARAMETR:  7.0794E-01  5.9092E-02  8.9762E-02  6.5904E-01  1.4316E+01  1.3433E+00  2.4107E-02  1.0000E-02  1.0000E-02  1.4638E+00
             1.0621E+01
 PARAMETER: -2.4539E-01 -2.7287E+00 -2.3106E+00 -3.1697E-01  2.7614E+00  3.9512E-01 -3.6252E+00 -1.9959E+01 -6.5428E+00  4.8106E-01
             2.4628E+00
 GRADIENT:  -7.6163E+00 -4.6105E-01 -2.3444E+01  4.9011E+01  2.0788E-01  1.0066E+01  7.5824E-04  0.0000E+00  0.0000E+00  1.9384E-02
            -2.0600E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -759.172703274925        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      696
 NPARAMETR:  6.5199E-01  2.8599E-02  7.0927E-02  5.5362E-01  1.1308E+01  1.2763E+00  1.0000E-02  1.0000E-02  1.0000E-02  1.0766E+00
             1.0813E+01
 PARAMETER: -3.2772E-01 -3.4544E+00 -2.5461E+00 -4.9128E-01  2.5255E+00  3.4396E-01 -4.8650E+00 -2.3316E+01 -8.3826E+00  1.7378E-01
             2.4807E+00
 GRADIENT:  -2.6472E+00  1.6139E-01 -2.2658E-02  1.2127E+00 -2.4046E-02 -7.0965E-01  0.0000E+00  0.0000E+00  0.0000E+00  1.2237E-04
             3.6910E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -759.189594415872        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      872
 NPARAMETR:  6.5112E-01  2.5188E-02  7.0313E-02  5.5049E-01  1.2943E+01  1.2788E+00  1.0000E-02  1.0000E-02  1.0000E-02  1.0670E+00
             1.0816E+01
 PARAMETER: -3.2907E-01 -3.5814E+00 -2.5548E+00 -4.9695E-01  2.6605E+00  3.4595E-01 -4.9068E+00 -2.3429E+01 -8.4293E+00  1.6490E-01
             2.4810E+00
 GRADIENT:   1.5188E+00  2.1793E-02 -4.8159E-01 -1.6421E-01  3.2585E-03 -2.8717E-02  0.0000E+00  0.0000E+00  0.0000E+00 -2.3102E-04
            -2.8524E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -759.191968789165        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     1058             RESET HESSIAN, TYPE I
 NPARAMETR:  6.4983E-01  2.2627E-02  7.0318E-02  5.5079E-01  1.2576E+01  1.2794E+00  1.0000E-02  1.0000E-02  1.0000E-02  1.1979E+00
             1.0811E+01
 PARAMETER: -3.3104E-01 -3.6886E+00 -2.5547E+00 -4.9640E-01  2.6318E+00  3.4640E-01 -4.9068E+00 -2.3429E+01 -8.4293E+00  2.8061E-01
             2.4805E+00
 GRADIENT:   2.1171E+01  3.0134E-02  2.1373E+01  1.3454E+01  7.7385E-03  3.0139E+00  0.0000E+00  0.0000E+00  0.0000E+00 -3.1386E-04
             1.9395E+01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -759.193057666059        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1237
 NPARAMETR:  6.5020E-01  2.0784E-02  7.0301E-02  5.5109E-01  9.5367E+00  1.2801E+00  1.0000E-02  1.0000E-02  1.0000E-02  1.4958E+00
             1.0820E+01
 PARAMETER: -3.3048E-01 -3.7736E+00 -2.5550E+00 -4.9586E-01  2.3551E+00  3.4697E-01 -4.9068E+00 -2.3429E+01 -8.4293E+00  5.0263E-01
             2.4814E+00
 GRADIENT:   1.0844E-02 -6.9993E-03  5.0454E-01 -1.0357E+00  1.1340E-03 -1.6359E-02  0.0000E+00  0.0000E+00  0.0000E+00 -5.2235E-04
             5.8037E-01

0ITERATION NO.:   53    OBJECTIVE VALUE:  -759.197068491932        NO. OF FUNC. EVALS.: 105
 CUMULATIVE NO. OF FUNC. EVALS.:     1342
 NPARAMETR:  6.5038E-01  2.0189E-02  7.0388E-02  5.5109E-01  9.3815E+00  1.2806E+00  1.0000E-02  1.0000E-02  1.0000E-02  2.6583E+00
             1.0819E+01
 PARAMETER: -3.3026E-01 -3.7927E+00 -2.5551E+00 -4.9470E-01  2.3156E+00  3.4764E-01 -4.9068E+00 -2.3429E+01 -8.4293E+00  1.0797E+00
             2.4803E+00
 GRADIENT:  -4.2687E-02  2.6595E-03 -5.9940E-01  1.3540E+00 -1.9432E+01  4.5033E-02  0.0000E+00  0.0000E+00  0.0000E+00  1.0580E-03
            -4.3113E-01

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1342
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.5604E-03  1.6722E-06  1.6902E-04 -3.3080E-04 -1.1802E-03
 SE:             2.8356E-02  2.0512E-06  2.5265E-04  3.9277E-04  1.0824E-03
 N:                     100         100         100         100         100

 P VAL.:         9.0008E-01  4.1494E-01  5.0349E-01  3.9966E-01  2.7556E-01

 ETASHRINKSD(%)  5.0048E+00  9.9993E+01  9.9154E+01  9.8684E+01  9.6374E+01
 ETASHRINKVR(%)  9.7590E+00  1.0000E+02  9.9993E+01  9.9983E+01  9.9869E+01
 EBVSHRINKSD(%)  5.1204E+00  9.9993E+01  9.9140E+01  9.8695E+01  9.7113E+01
 EBVSHRINKVR(%)  9.9787E+00  1.0000E+02  9.9993E+01  9.9983E+01  9.9917E+01
 RELATIVEINF(%)  5.7051E+00  1.1985E-08  4.8625E-05  7.0181E-05  9.1991E-03
 EPSSHRINKSD(%)  6.5315E+00
 EPSSHRINKVR(%)  1.2636E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -759.19706849193187     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -24.046241928193695     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    17.46
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.93
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -759.197       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         6.50E-01  2.04E-02  7.03E-02  5.52E-01  9.17E+00  1.28E+00  1.00E-02  1.00E-02  1.00E-02  2.66E+00  1.08E+01
 


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
+        1.42E+03
 
 TH 2
+       -4.16E+02  5.80E+04
 
 TH 3
+       -1.83E+03 -2.46E+04  1.02E+05
 
 TH 4
+       -4.51E+02  3.49E+02 -2.09E+04  3.89E+03
 
 TH 5
+        1.25E+01 -2.93E+02 -2.59E+00  5.09E+01  3.95E+00
 
 TH 6
+       -6.27E-01  1.05E+01  2.24E+02 -5.90E+01 -6.61E-01  9.49E+01
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH10
+        4.56E-02  6.54E+00 -1.25E+02  2.08E-01  2.25E+00  5.27E-03  0.00E+00  0.00E+00  0.00E+00  8.86E-02
 
 TH11
+       -2.04E+01 -1.18E+02  1.27E+02 -2.11E+01 -2.58E-01  1.14E+00  0.00E+00  0.00E+00  0.00E+00 -1.30E-02  4.12E+00
 
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
 #CPUT: Total CPU Time in Seconds,       24.447
Stop Time:
Wed Sep 29 20:02:30 CDT 2021
