Sat Sep 25 02:19:37 CDT 2021
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
$DATA ../../../../data/int/SL3/dat47.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      981
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

 TOT. NO. OF OBS RECS:      881
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
 RAW OUTPUT FILE (FILE): m47.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -430.090905658038        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.4788E+01 -1.3055E+02  1.7139E+02  6.4603E+01  1.7871E+02  3.4484E+01 -1.5004E+02 -2.4708E+02 -1.6345E+02 -5.5999E+01
            -6.1943E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2682.26445618003        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.0240E+00  1.3547E+00  8.5262E-01  8.7848E-01  1.0542E+00  8.1363E-01  1.1202E+00  1.0181E+00  1.1736E+00  1.0889E+00
             2.7956E+00
 PARAMETER:  1.2374E-01  4.0357E-01 -5.9439E-02 -2.9557E-02  1.5279E-01 -1.0625E-01  2.1354E-01  1.1794E-01  2.6009E-01  1.8516E-01
             1.1280E+00
 GRADIENT:   4.0032E+01  6.0439E+00 -1.0610E+01  1.3923E+01 -6.7205E+00 -4.0581E+01  2.0727E+01  3.6083E+00 -1.9414E+00 -1.0980E+01
             6.6333E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2689.51318985588        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  1.0272E+00  1.8140E+00  8.7959E-01  7.1048E-01  1.2833E+00  8.7918E-01  8.4712E-01  9.3564E-01  1.6208E+00  1.4501E+00
             2.7649E+00
 PARAMETER:  1.2681E-01  6.9551E-01 -2.8301E-02 -2.4182E-01  3.4940E-01 -2.8771E-02 -6.5909E-02  3.3477E-02  5.8290E-01  4.7161E-01
             1.1170E+00
 GRADIENT:   4.2855E+01  1.4893E+02  5.6796E+00  8.5819E+01 -5.0011E+01 -1.0425E+01  7.7912E+00  6.5113E-01  1.0649E+01  2.6912E+00
             5.5723E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2707.09256956352        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      228
 NPARAMETR:  1.0047E+00  2.0874E+00  8.1446E-01  4.4829E-01  1.6069E+00  9.0694E-01  6.5822E-01  5.5038E-01  1.8916E+00  1.6643E+00
             2.6438E+00
 PARAMETER:  1.0470E-01  8.3593E-01 -1.0523E-01 -7.0231E-01  5.7432E-01  2.3260E-03 -3.1821E-01 -4.9714E-01  7.3743E-01  6.0938E-01
             1.0722E+00
 GRADIENT:  -1.2235E+01  1.0474E+02  3.4533E+00  3.4489E+01  3.0930E+00  9.7053E-01 -1.1608E+01 -3.1536E-01 -1.2784E+01  6.4264E+00
            -2.6371E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2716.15785946630        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      300
 NPARAMETR:  1.0101E+00  2.3307E+00  4.7660E-01  2.3449E-01  1.7605E+00  8.9494E-01  6.6982E-01  6.2492E-02  2.9506E+00  1.7400E+00
             2.6337E+00
 PARAMETER:  1.1006E-01  9.4618E-01 -6.4108E-01 -1.3503E+00  6.6558E-01 -1.0999E-02 -3.0075E-01 -2.6727E+00  1.1820E+00  6.5388E-01
             1.0684E+00
 GRADIENT:   4.9916E+00  2.5099E+01 -1.4392E+00  8.5071E+00  4.7197E+00 -3.6197E+00  2.1975E+00  4.8117E-03  3.2973E+00  2.6070E-01
            -8.3844E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2718.62573266486        NO. OF FUNC. EVALS.:  81
 CUMULATIVE NO. OF FUNC. EVALS.:      381
 NPARAMETR:  1.0031E+00  2.4940E+00  3.8175E-01  1.1607E-01  1.8662E+00  9.1863E-01  6.2811E-01  1.0000E-02  4.1621E+00  1.8280E+00
             2.6400E+00
 PARAMETER:  1.0305E-01  1.0139E+00 -8.6299E-01 -2.0536E+00  7.2393E-01  1.5129E-02 -3.6503E-01 -4.6374E+00  1.5260E+00  7.0324E-01
             1.0708E+00
 GRADIENT:  -1.3147E+01  4.2900E+01 -5.0971E-01  1.4111E+00  2.3020E+00  5.8775E+00 -8.1422E+00  0.0000E+00 -1.0552E+00  1.0942E+00
             7.9531E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2718.93441481118        NO. OF FUNC. EVALS.: 115
 CUMULATIVE NO. OF FUNC. EVALS.:      496
 NPARAMETR:  1.0080E+00  2.5000E+00  3.7135E-01  1.1250E-01  1.8603E+00  9.0957E-01  6.4785E-01  1.0000E-02  4.2555E+00  1.8325E+00
             2.6348E+00
 PARAMETER:  1.0796E-01  1.0163E+00 -8.9061E-01 -2.0848E+00  7.2072E-01  5.2143E-03 -3.3410E-01 -4.7054E+00  1.5482E+00  7.0567E-01
             1.0688E+00
 GRADIENT:  -6.7623E+00  1.9396E+00  5.4367E-01 -3.8203E-01 -3.8290E+00  1.8184E+00 -1.2109E+00  0.0000E+00 -1.8192E+00  1.0302E-01
             2.4427E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2719.01253375847        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      671
 NPARAMETR:  1.0106E+00  2.5126E+00  3.4893E-01  1.0448E-01  1.8811E+00  9.0517E-01  6.5172E-01  1.0000E-02  4.4504E+00  1.8425E+00
             2.6310E+00
 PARAMETER:  1.1052E-01  1.0213E+00 -9.5288E-01 -2.1588E+00  7.3185E-01  3.6514E-04 -3.2814E-01 -5.0058E+00  1.5930E+00  7.1114E-01
             1.0674E+00
 GRADIENT:  -1.0722E-02  8.2343E-02  3.2943E-02 -3.0970E-02  3.8733E-02 -3.6859E-03  1.2476E-02  0.0000E+00 -2.1345E-02  9.4994E-03
            -2.5736E-02

0ITERATION NO.:   36    OBJECTIVE VALUE:  -2719.01253375847        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:      693
 NPARAMETR:  1.0106E+00  2.5126E+00  3.4893E-01  1.0448E-01  1.8811E+00  9.0517E-01  6.5172E-01  1.0000E-02  4.4504E+00  1.8425E+00
             2.6310E+00
 PARAMETER:  1.1052E-01  1.0213E+00 -9.5288E-01 -2.1588E+00  7.3185E-01  3.6514E-04 -3.2814E-01 -5.0058E+00  1.5930E+00  7.1114E-01
             1.0674E+00
 GRADIENT:  -1.0722E-02  8.2343E-02  3.2943E-02 -3.0970E-02  3.8733E-02 -3.6859E-03  1.2476E-02  0.0000E+00 -2.1345E-02  9.4994E-03
            -2.5736E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      693
 NO. OF SIG. DIGITS IN FINAL EST.:  3.6
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.9259E-03 -4.2713E-02  1.3452E-05  5.1342E-02 -2.7775E-02
 SE:             2.9354E-02  2.4014E-02  2.3410E-05  1.8084E-02  2.6056E-02
 N:                     100         100         100         100         100

 P VAL.:         9.4769E-01  7.5296E-02  5.6554E-01  4.5248E-03  2.8643E-01

 ETASHRINKSD(%)  1.6589E+00  1.9549E+01  9.9922E+01  3.9416E+01  1.2709E+01
 ETASHRINKVR(%)  3.2903E+00  3.5276E+01  1.0000E+02  6.3296E+01  2.3803E+01
 EBVSHRINKSD(%)  1.9271E+00  1.4579E+01  9.9886E+01  4.7899E+01  8.5792E+00
 EBVSHRINKVR(%)  3.8170E+00  2.7032E+01  1.0000E+02  7.2855E+01  1.6422E+01
 RELATIVEINF(%)  9.6112E+01  3.0907E+01  1.0321E-04  1.1776E+01  7.0159E+01
 EPSSHRINKSD(%)  1.6466E+01
 EPSSHRINKVR(%)  3.0221E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          881
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1619.1696955066332     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2719.0125337584745     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1099.8428382518414     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    14.07
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    13.33
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2719.013       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.01E+00  2.51E+00  3.49E-01  1.04E-01  1.88E+00  9.05E-01  6.52E-01  1.00E-02  4.45E+00  1.84E+00  2.63E+00
 


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
+        1.27E+03
 
 TH 2
+       -1.97E+01  3.68E+02
 
 TH 3
+       -1.13E+01  1.31E+02  2.95E+02
 
 TH 4
+        1.06E+01 -2.35E+01 -8.20E+02  4.50E+03
 
 TH 5
+       -3.30E+00 -1.53E+01  8.11E+00  3.94E+00  7.80E+01
 
 TH 6
+        3.04E+00 -5.49E+00 -1.84E+00  9.98E+00 -1.19E+00  2.19E+02
 
 TH 7
+        3.88E+00 -6.59E-02  8.16E+01 -2.44E+02 -6.16E+00 -3.24E+00  2.95E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        9.02E-01 -9.30E+00 -1.37E+01  9.68E+01 -5.94E-01  1.57E-01 -6.89E-01  0.00E+00  3.92E+00
 
 TH10
+        3.72E-01 -2.85E+00  4.99E+00  1.52E+01 -6.08E+00  5.99E-01 -4.69E+00  0.00E+00  2.16E-02  3.89E+01
 
 TH11
+       -1.77E+01 -1.39E+01  2.02E-01 -2.55E+01  1.05E+00  2.81E+00  7.88E+00  0.00E+00  3.42E-01  4.95E+00  1.66E+02
 
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
 #CPUT: Total CPU Time in Seconds,       27.508
Stop Time:
Sat Sep 25 02:20:06 CDT 2021
