Sat Sep 18 11:20:05 CDT 2021
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
$DATA ../../../../data/spa/S1/dat85.csv ignore=@
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
Current Date:       18 SEP 2021
Days until program expires : 211
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
 (E4.0,E3.0,E20.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m85.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1179.34528984827        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   8.5055E+01 -8.4646E+01  4.2005E+01 -9.9399E+01  2.4740E+01  2.3260E+01 -2.1964E+01 -1.9581E+02 -2.4227E+01  1.7044E+01
            -7.2178E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1655.11241702301        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:       89
 NPARAMETR:  1.0269E+00  1.1456E+00  1.1405E+00  9.8348E-01  1.0399E+00  9.0061E-01  1.1574E+00  1.0339E+00  1.0531E+00  8.6061E-01
             8.8637E-01
 PARAMETER:  1.2658E-01  2.3590E-01  2.3150E-01  8.3345E-02  1.3913E-01 -4.6828E-03  2.4617E-01  1.3330E-01  1.5170E-01 -5.0114E-02
            -2.0624E-02
 GRADIENT:   1.9325E+02  5.0186E+01  5.3172E+01  6.3602E-01 -7.9918E+01 -2.7218E+01 -8.2647E+00 -1.0228E+01 -1.0322E+00 -4.1683E+00
            -5.6210E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1659.22990926916        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      162
 NPARAMETR:  1.0180E+00  9.7737E-01  9.4993E-01  1.0879E+00  8.9772E-01  9.4278E-01  1.4498E+00  9.6796E-01  9.5605E-01  5.7377E-01
             9.0452E-01
 PARAMETER:  1.1781E-01  7.7112E-02  4.8632E-02  1.8423E-01 -7.8951E-03  4.1081E-02  4.7145E-01  6.7431E-02  5.5052E-02 -4.5553E-01
            -3.4805E-04
 GRADIENT:   1.5089E+02  3.6218E+01  3.4691E+01  4.1307E+01 -4.5649E+01 -5.6371E+00 -3.3886E+00 -1.4870E+00  6.6493E+00 -1.3175E+01
            -4.6529E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1666.76469971123        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      233
 NPARAMETR:  9.7441E-01  9.6930E-01  7.7320E-01  1.0486E+00  8.4480E-01  9.4113E-01  1.4861E+00  5.8132E-01  9.1048E-01  6.5283E-01
             9.7417E-01
 PARAMETER:  7.4078E-02  6.8822E-02 -1.5722E-01  1.4741E-01 -6.8656E-02  3.9325E-02  4.9615E-01 -4.4246E-01  6.2177E-03 -3.2644E-01
             7.3833E-02
 GRADIENT:   2.2249E+01  2.2998E-01 -4.5898E+00  1.2764E+01  2.9540E+00 -3.3962E+00  2.4594E+00  2.0253E+00  3.2503E+00  1.0814E+00
            -5.6746E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1667.27151479236        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:      394
 NPARAMETR:  9.8089E-01  1.0391E+00  7.0781E-01  1.0044E+00  8.3749E-01  9.5416E-01  1.4268E+00  3.8459E-01  9.2725E-01  6.3931E-01
             9.8770E-01
 PARAMETER:  8.0709E-02  1.3836E-01 -2.4559E-01  1.0440E-01 -7.7350E-02  5.3076E-02  4.5542E-01 -8.5558E-01  2.4466E-02 -3.4737E-01
             8.7619E-02
 GRADIENT:  -4.5616E+00  2.5102E+00 -2.6459E+00  4.9845E+00  2.0880E+00 -3.1193E+00  1.8224E+00  5.7255E-01  1.2306E+00  3.1654E-01
            -2.9280E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1667.46588145237        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      570
 NPARAMETR:  9.8311E-01  1.1116E+00  6.5305E-01  9.5144E-01  8.3975E-01  9.6096E-01  1.3345E+00  1.6812E-01  9.5714E-01  6.3091E-01
             9.8897E-01
 PARAMETER:  8.2970E-02  2.0577E-01 -3.2610E-01  5.0220E-02 -7.4649E-02  6.0180E-02  3.8857E-01 -1.6831E+00  5.6190E-02 -3.6060E-01
             8.8909E-02
 GRADIENT:  -4.6270E-01 -6.9630E-01 -6.9319E-01 -9.1172E-01  4.2703E-01 -5.3991E-01 -3.6619E-02  1.1129E-01  4.4166E-01  1.4329E-01
             3.4524E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1667.51998332354        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      746
 NPARAMETR:  9.8371E-01  1.0836E+00  6.4640E-01  9.6730E-01  8.2201E-01  9.6234E-01  1.3627E+00  2.3694E-02  9.3897E-01  6.2126E-01
             9.8863E-01
 PARAMETER:  8.3573E-02  1.8029E-01 -3.3634E-01  6.6754E-02 -9.6002E-02  6.1615E-02  4.0949E-01 -3.6425E+00  3.7027E-02 -3.7601E-01
             8.8564E-02
 GRADIENT:   9.4807E-01 -3.1411E-02 -3.6641E-01  4.3007E-01  1.8918E-01  2.4574E-02 -2.5726E-01  1.8218E-03 -1.9336E-01  6.8158E-02
             1.1172E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1667.52154187392        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      922
 NPARAMETR:  9.8331E-01  1.0848E+00  6.4783E-01  9.6654E-01  8.2357E-01  9.6226E-01  1.3634E+00  1.0000E-02  9.4098E-01  6.2238E-01
             9.8843E-01
 PARAMETER:  8.3172E-02  1.8141E-01 -3.3413E-01  6.5970E-02 -9.4105E-02  6.1534E-02  4.1001E-01 -4.7279E+00  3.9166E-02 -3.7420E-01
             8.8358E-02
 GRADIENT:   3.6520E-02 -4.5200E-03 -6.7595E-03  4.6512E-03  9.9657E-03 -3.2820E-03 -5.4833E-03  0.0000E+00 -5.4528E-03  3.6099E-04
             5.8493E-03

0ITERATION NO.:   36    OBJECTIVE VALUE:  -1667.52154187392        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:      944
 NPARAMETR:  9.8331E-01  1.0848E+00  6.4783E-01  9.6654E-01  8.2357E-01  9.6226E-01  1.3634E+00  1.0000E-02  9.4098E-01  6.2238E-01
             9.8843E-01
 PARAMETER:  8.3172E-02  1.8141E-01 -3.3413E-01  6.5970E-02 -9.4105E-02  6.1534E-02  4.1001E-01 -4.7279E+00  3.9166E-02 -3.7420E-01
             8.8358E-02
 GRADIENT:   3.6520E-02 -4.5200E-03 -6.7595E-03  4.6512E-03  9.9657E-03 -3.2820E-03 -5.4833E-03  0.0000E+00 -5.4528E-03  3.6099E-04
             5.8493E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      944
 NO. OF SIG. DIGITS IN FINAL EST.:  3.3
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.2725E-04  1.8037E-03 -5.0830E-04 -2.9898E-03 -8.3920E-03
 SE:             2.9849E-02  2.4709E-02  2.0179E-04  2.4534E-02  1.9167E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9660E-01  9.4181E-01  1.1769E-02  9.0301E-01  6.6151E-01

 ETASHRINKSD(%)  3.6647E-03  1.7222E+01  9.9324E+01  1.7810E+01  3.5788E+01
 ETASHRINKVR(%)  7.3292E-03  3.1478E+01  9.9995E+01  3.2447E+01  5.8768E+01
 EBVSHRINKSD(%)  4.4309E-01  1.6132E+01  9.9396E+01  1.8051E+01  3.6523E+01
 EBVSHRINKVR(%)  8.8422E-01  2.9661E+01  9.9996E+01  3.2844E+01  5.9707E+01
 RELATIVEINF(%)  9.8990E+01  7.9288E+00  4.2231E-04  8.0385E+00  3.5115E+00
 EPSSHRINKSD(%)  4.3654E+01
 EPSSHRINKVR(%)  6.8251E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1667.5215418739178     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -932.37071531017966     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    10.34
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.55
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1667.522       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.83E-01  1.08E+00  6.48E-01  9.67E-01  8.24E-01  9.62E-01  1.36E+00  1.00E-02  9.41E-01  6.22E-01  9.88E-01
 


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
+        1.23E+03
 
 TH 2
+       -4.66E+00  3.67E+02
 
 TH 3
+        1.68E+01  2.85E+02  1.07E+03
 
 TH 4
+       -1.12E+01  2.02E+02 -5.13E+02  9.20E+02
 
 TH 5
+       -7.61E+00 -4.27E+02 -1.18E+03  5.95E+02  1.71E+03
 
 TH 6
+        1.59E+00 -9.06E-01  2.85E+00 -3.52E+00 -1.93E+00  2.12E+02
 
 TH 7
+        8.05E-01  2.79E+01 -4.08E+01 -1.55E+01  2.30E+01  1.27E-01  5.53E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        2.09E+00 -1.98E+01 -3.97E+01  4.52E+01 -2.09E+01 -4.50E-01  8.85E+00  0.00E+00  1.14E+02
 
 TH10
+       -1.46E+00 -1.37E+01 -7.80E+01 -3.56E+01 -6.45E+01 -1.55E-01  1.63E+01  0.00E+00  1.93E+01  1.04E+02
 
 TH11
+       -6.64E+00 -1.16E+01 -3.81E+01 -6.59E+00 -2.30E+00  1.67E+00  5.69E+00  0.00E+00  1.42E+01  2.69E+01  2.21E+02
 
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
 #CPUT: Total CPU Time in Seconds,       15.952
Stop Time:
Sat Sep 18 11:20:23 CDT 2021
