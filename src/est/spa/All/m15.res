Sat Sep 25 14:55:52 CDT 2021
$PROB template control stream
;-----------------------------------------------------------------------
; Project: 	Investigating the contribution of residual unexplained
; 	   	variability in nonlinear mixed-effect approach
; Model: 	One-compartment model with linear elimination
; Estim:	First-order conditional est. with interaction
; Author: 	Mutaz M. Jaber <jaber038@umn.edu>
; Date created: 9/7/2021
; Date modified: 9/7/2021
;-----------------------------------------------------------------------
$INPUT ID TIME DV AMT MDV EVID
$DATA ../../../../data/spa/All/dat15.csv ignore=@
$SUBR ADVAN2 TRANS2
$EST MET=1 NOABORT MAX=10000 PRINT=5 INTER
$PK

ET1 = EXP(ETA(1)*THETA(4))
ET2 = EXP(ETA(2)*THETA(5))
ET3 = EXP(ETA(3)*THETA(6))


CL = 5.0 * THETA(1) * ET1
V = 85  * THETA(2) * ET2
KA = 0.7 * THETA(3) * ET3

SC = V
$ERROR
CVERR 	= 0.05
W  	= THETA(7)*F*CVERR
Y  	= F + W * ERR(1)
$THETA
(0,1) ; tvCL
(0,1) ; tvV
(0,1) ; tvKA
(0,1) ; tvCL
(0,1) ; tvV
(0,1) ; tvK
(0,1) ; RUV
$OMEGA
0.9 FIX ;     IIV CL
0.9 FIX  ;     IIV V
0.9 FIX ;      IIV KA
$SIGMA  1  FIX;        [P]
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
0LENGTH OF THETA:   7
0DEFAULT THETA BOUNDARY TEST OMITTED:    NO
0OMEGA HAS SIMPLE DIAGONAL FORM WITH DIMENSION:   3
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
0INITIAL ESTIMATE OF OMEGA:
 0.9000E+00
 0.0000E+00   0.9000E+00
 0.0000E+00   0.0000E+00   0.9000E+00
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

 ONE COMPARTMENT MODEL WITH FIRST-ORDER ABSORPTION (ADVAN2)
0MAXIMUM NO. OF BASIC PK PARAMETERS:   3
0BASIC PK PARAMETERS (AFTER TRANSLATION):
   ELIMINATION RATE (K) IS BASIC PK PARAMETER NO.:  1
   ABSORPTION RATE (KA) IS BASIC PK PARAMETER NO.:  3

 TRANSLATOR WILL CONVERT PARAMETERS
 CLEARANCE (CL) AND VOLUME (V) TO K (TRANS2)
0COMPARTMENT ATTRIBUTES
 COMPT. NO.   FUNCTION   INITIAL    ON/OFF      DOSE      DEFAULT    DEFAULT
                         STATUS     ALLOWED    ALLOWED    FOR DOSE   FOR OBS.
    1         DEPOT        OFF        YES        YES        YES        NO
    2         CENTRAL      ON         NO         YES        NO         YES
    3         OUTPUT       OFF        YES        NO         NO         NO
1
 ADDITIONAL PK PARAMETERS - ASSIGNMENT OF ROWS IN GG
 COMPT. NO.                             INDICES
              SCALE      BIOAVAIL.   ZERO-ORDER  ZERO-ORDER  ABSORB
                         FRACTION    RATE        DURATION    LAG
    1            *           *           *           *           *
    2            4           *           *           *           *
    3            *           -           -           -           -
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
 RAW OUTPUT FILE (FILE): m15.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   808.648187102694        NO. OF FUNC. EVALS.:   9
 CUMULATIVE NO. OF FUNC. EVALS.:        9
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
 GRADIENT:   3.2284E+01  1.0165E+01 -7.4390E+01  1.5180E+02  1.2298E+02 -7.5118E+00 -4.1759E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -876.294090110689        NO. OF FUNC. EVALS.:  51
 CUMULATIVE NO. OF FUNC. EVALS.:       60
 NPARAMETR:  1.3109E+00  1.0742E+00  2.2305E+00  3.4489E-01  2.5040E-01  5.8087E-01  1.0288E+01
 PARAMETER:  3.7071E-01  1.7154E-01  9.0222E-01 -9.6452E-01 -1.2847E+00 -4.4323E-01  2.4310E+00
 GRADIENT:   1.6264E+02 -2.4476E+02  3.9474E+00 -4.3932E+00  8.3602E+00  5.6405E+00  2.5687E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -971.362179958778        NO. OF FUNC. EVALS.:  52
 CUMULATIVE NO. OF FUNC. EVALS.:      112
 NPARAMETR:  1.0951E+00  1.2187E+00  4.4123E+01  4.3338E-01  3.5566E-01  2.5476E+00  6.6731E+00
 PARAMETER:  1.9082E-01  2.9778E-01  3.8870E+00 -7.3615E-01 -9.3377E-01  1.0352E+00  1.9981E+00
 GRADIENT:  -5.2507E+01 -1.8117E+01 -1.4527E-01  8.6931E+00 -6.3985E+00 -1.8197E-02  8.7154E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -984.679917532348        NO. OF FUNC. EVALS.:  50
 CUMULATIVE NO. OF FUNC. EVALS.:      162
 NPARAMETR:  1.1477E+00  1.2610E+00  1.3277E+02  4.8013E-01  4.8555E-01  2.7460E+00  5.4076E+00
 PARAMETER:  2.3776E-01  3.3189E-01  4.9886E+00 -6.3370E-01 -6.2247E-01  1.1102E+00  1.7878E+00
 GRADIENT:   4.5131E+00 -2.6842E+00  4.5911E-03  9.0869E-01  6.2293E-01 -1.5509E-03  1.0742E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -984.680609731897        NO. OF FUNC. EVALS.:  63
 CUMULATIVE NO. OF FUNC. EVALS.:      225
 NPARAMETR:  1.1443E+00  1.2625E+00  1.3435E+02  4.7919E-01  4.8463E-01  2.7656E+00  5.4044E+00
 PARAMETER:  2.3481E-01  3.3306E-01  5.0005E+00 -6.3565E-01 -6.2436E-01  1.1173E+00  1.7872E+00
 GRADIENT:  -1.5186E+00 -3.4198E+00  3.5735E-03 -1.4640E+00 -5.7237E-01 -1.5775E-03 -9.9083E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -984.700606208395        NO. OF FUNC. EVALS.: 108
 CUMULATIVE NO. OF FUNC. EVALS.:      333
 NPARAMETR:  1.1465E+00  1.2689E+00  1.3575E+02  4.8165E-01  4.8646E-01  8.2704E+00  5.4075E+00
 PARAMETER:  2.3671E-01  3.3815E-01  5.0108E+00 -6.3054E-01 -6.2059E-01  2.2127E+00  1.7878E+00
 GRADIENT:  -3.5268E-01 -4.8666E-02  9.7547E-03 -1.2385E-01 -5.3630E-02 -5.6346E-04 -1.9621E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -984.712087659359        NO. OF FUNC. EVALS.: 125
 CUMULATIVE NO. OF FUNC. EVALS.:      458
 NPARAMETR:  1.1525E+00  1.2720E+00  3.5785E+01  4.8236E-01  4.8583E-01  4.8099E+00  5.4243E+00
 PARAMETER:  2.4191E-01  3.4062E-01  3.6775E+00 -6.2905E-01 -6.2189E-01  1.6707E+00  1.7909E+00
 GRADIENT:   2.6441E+00  2.4522E-01 -2.3654E-01  4.4883E-01 -2.4088E-01  5.4087E-01  1.6428E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -984.756327578487        NO. OF FUNC. EVALS.: 120
 CUMULATIVE NO. OF FUNC. EVALS.:      578
 NPARAMETR:  1.1576E+00  1.2751E+00  1.2553E+01  4.8270E-01  4.8486E-01  2.8294E+00  5.4385E+00
 PARAMETER:  2.4636E-01  3.4301E-01  2.6300E+00 -6.2836E-01 -6.2390E-01  1.1401E+00  1.7935E+00
 GRADIENT:   5.0995E+00  7.8408E-02 -3.1694E+00  8.5231E-01 -2.2019E+00  5.9606E+00  1.3058E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -984.835687235456        NO. OF FUNC. EVALS.: 119
 CUMULATIVE NO. OF FUNC. EVALS.:      697
 NPARAMETR:  1.1563E+00  1.2801E+00  8.4925E+00  4.8217E-01  4.8439E-01  2.2315E+00  5.4109E+00
 PARAMETER:  2.4521E-01  3.4696E-01  2.2392E+00 -6.2947E-01 -6.2487E-01  9.0268E-01  1.7884E+00
 GRADIENT:   1.1831E-01  9.5641E-02 -1.7958E+00 -7.8093E-02 -7.3582E-01  3.0455E+00 -1.2996E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -984.876145972229        NO. OF FUNC. EVALS.: 116
 CUMULATIVE NO. OF FUNC. EVALS.:      813
 NPARAMETR:  1.1601E+00  1.2842E+00  6.9895E+00  4.8327E-01  4.8492E-01  1.9492E+00  5.3956E+00
 PARAMETER:  2.4847E-01  3.5010E-01  2.0444E+00 -6.2718E-01 -6.2377E-01  7.6740E-01  1.7856E+00
 GRADIENT:   9.9198E-01  5.9809E-01 -1.7375E+00  1.1140E-01 -3.7515E-01  2.6642E+00 -3.1330E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -984.899484440958        NO. OF FUNC. EVALS.: 116
 CUMULATIVE NO. OF FUNC. EVALS.:      929
 NPARAMETR:  1.1622E+00  1.2835E+00  6.6543E+00  4.8344E-01  4.8614E-01  1.8591E+00  5.4174E+00
 PARAMETER:  2.5030E-01  3.4959E-01  1.9953E+00 -6.2682E-01 -6.2126E-01  7.2007E-01  1.7896E+00
 GRADIENT:   1.7104E+00 -5.4909E-01 -1.6699E+00  1.5740E-01  5.6334E-01  1.9825E+00  8.4546E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1001.11613747105        NO. OF FUNC. EVALS.: 128
 CUMULATIVE NO. OF FUNC. EVALS.:     1057
 NPARAMETR:  1.1582E+00  1.2255E+00  1.4005E+00  4.7692E-01  5.1574E-01  1.9821E-01  5.1215E+00
 PARAMETER:  2.4689E-01  3.0336E-01  4.3683E-01 -6.4041E-01 -5.6214E-01 -1.5184E+00  1.7334E+00
 GRADIENT:  -1.7850E+00  6.3860E+00 -3.2710E+00 -1.1378E+01  7.6863E+00  2.6896E+00  8.3475E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1001.65352673906        NO. OF FUNC. EVALS.: 116
 CUMULATIVE NO. OF FUNC. EVALS.:     1173
 NPARAMETR:  1.1463E+00  1.2073E+00  1.3821E+00  4.9181E-01  5.0679E-01  1.7334E-01  5.0037E+00
 PARAMETER:  2.3657E-01  2.8841E-01  4.2359E-01 -6.0966E-01 -5.7965E-01 -1.6525E+00  1.7102E+00
 GRADIENT:  -7.6930E+00 -1.0348E+00 -3.1361E+00 -1.0928E+00  1.5354E-01  2.0125E+00 -7.2850E+00

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1002.67385703785        NO. OF FUNC. EVALS.: 117
 CUMULATIVE NO. OF FUNC. EVALS.:     1290
 NPARAMETR:  1.1587E+00  1.2125E+00  1.4000E+00  4.9195E-01  5.0401E-01  3.8427E-02  5.0720E+00
 PARAMETER:  2.4726E-01  2.9269E-01  4.3648E-01 -6.0937E-01 -5.8517E-01 -3.1590E+00  1.7237E+00
 GRADIENT:   4.0233E-01  4.0879E-02 -3.9005E-01 -2.6796E-01 -1.8978E-01  8.8618E-02 -1.0686E-01

0ITERATION NO.:   69    OBJECTIVE VALUE:  -1002.71560593506        NO. OF FUNC. EVALS.:  83
 CUMULATIVE NO. OF FUNC. EVALS.:     1373
 NPARAMETR:  1.1584E+00  1.2129E+00  1.4023E+00  4.9234E-01  5.0434E-01  1.0000E-02  5.0735E+00
 PARAMETER:  2.4703E-01  2.9300E-01  4.3812E-01 -6.0858E-01 -5.8450E-01 -4.6308E+00  1.7240E+00
 GRADIENT:   1.3683E-01  7.8355E-02 -4.2601E-02 -1.8549E-02 -2.0460E-02  0.0000E+00  1.3015E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1373
 NO. OF SIG. DIGITS IN FINAL EST.:  3.2
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -2.6281E-03 -3.7869E-02 -7.6288E-04
 SE:             9.0977E-02  8.7845E-02  7.9966E-04
 N:                     100         100         100

 P VAL.:         9.7695E-01  6.6640E-01  3.4008E-01

 ETASHRINKSD(%)  3.6186E+00  6.9371E+00  9.9153E+01
 ETASHRINKVR(%)  7.1062E+00  1.3393E+01  9.9993E+01
 EBVSHRINKSD(%)  3.3321E+00  6.6089E+00  9.8988E+01
 EBVSHRINKVR(%)  6.5531E+00  1.2781E+01  9.9990E+01
 RELATIVEINF(%)  9.2439E+01  8.0241E+01  9.4329E-03
 EPSSHRINKSD(%)  2.1974E+01
 EPSSHRINKVR(%)  3.9120E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1002.7156059350555     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -267.56477937131729     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           300
  
 #TERE:
 Elapsed estimation  time in seconds:    12.77
0R MATRIX ALGORITHMICALLY SINGULAR
0COVARIANCE MATRIX UNOBTAINABLE
0R MATRIX IS OUTPUT
0S MATRIX ALGORITHMICALLY SINGULAR
0S MATRIX IS OUTPUT
0T MATRIX - EQUAL TO RS*RMAT, WHERE S* IS A PSEUDO INVERSE OF S - IS OUTPUT
 Elapsed covariance  time in seconds:     2.03
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1002.716       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7     
 
         1.16E+00  1.21E+00  1.40E+00  4.92E-01  5.04E-01  1.00E-02  5.07E+00
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


         ETA1      ETA2      ETA3     
 
 ETA1
+        9.00E-01
 
 ETA2
+        0.00E+00  9.00E-01
 
 ETA3
+        0.00E+00  0.00E+00  9.00E-01
 


 SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS  ****


         EPS1     
 
 EPS1
+        1.00E+00
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


         ETA1      ETA2      ETA3     
 
 ETA1
+        9.49E-01
 
 ETA2
+        0.00E+00  9.49E-01
 
 ETA3
+        0.00E+00  0.00E+00  9.49E-01
 


 SIGMA - CORR MATRIX FOR RANDOM EFFECTS - EPSILONS  ***


         EPS1     
 
 EPS1
+        1.00E+00
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                                     T MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      OM11      OM12      OM13      OM22      OM23  
             OM33      SG11  
 
 TH 1
+        4.37E+02
 
 TH 2
+       -2.28E+02  3.75E+02
 
 TH 3
+       -1.08E+01 -2.45E+01  3.89E+00
 
 TH 4
+        1.50E+02 -6.68E+01 -1.23E+01  1.02E+03
 
 TH 5
+       -1.08E+02 -1.08E+02  2.80E+01 -5.88E+02  7.19E+02
 
 TH 6
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 7
+        3.97E+00 -1.78E+01  1.85E+00 -1.52E+00  2.15E+01  0.00E+00  1.45E+00
 
 OM11
+       ......... ......... ......... ......... ......... ......... ......... .........
 
 OM12
+       ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM13
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM22
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM23
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM33
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         .........
 
 SG11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... .........
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                                     R MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      OM11      OM12      OM13      OM22      OM23  
             OM33      SG11  
 
 TH 1
+        3.18E+02
 
 TH 2
+       -2.12E+01  2.52E+02
 
 TH 3
+       -7.56E+00 -3.52E+01  5.07E+01
 
 TH 4
+       -9.69E+00 -2.05E+01  1.22E+00  6.97E+02
 
 TH 5
+       -2.51E+00 -1.66E+01  5.14E+00 -7.64E+01  5.16E+02
 
 TH 6
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 7
+       -6.62E+00 -5.02E+00 -9.47E-01  1.16E+01  1.83E+01  0.00E+00  1.48E+01
 
 OM11
+       ......... ......... ......... ......... ......... ......... ......... .........
 
 OM12
+       ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM13
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM22
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM23
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM33
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         .........
 
 SG11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... .........
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                                     S MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      OM11      OM12      OM13      OM22      OM23  
             OM33      SG11  
 
 TH 1
+        3.19E+02
 
 TH 2
+        1.98E+02  2.60E+02
 
 TH 3
+       -5.90E+01 -2.62E+01  4.15E+01
 
 TH 4
+        1.12E-01  7.89E+01  1.24E+00  8.30E+02
 
 TH 5
+        1.55E+02  1.87E+02 -1.58E+01  4.39E+02  7.60E+02
 
 TH 6
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 7
+        2.92E+01 -7.95E+00 -1.47E+01  1.60E+01  4.95E+01  0.00E+00  3.76E+01
 
 OM11
+       ......... ......... ......... ......... ......... ......... ......... .........
 
 OM12
+       ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM13
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM22
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM23
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM33
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         .........
 
 SG11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... .........
 
 Elapsed finaloutput time in seconds:     0.02
 #CPUT: Total CPU Time in Seconds,       14.836
Stop Time:
Sat Sep 25 14:56:11 CDT 2021
