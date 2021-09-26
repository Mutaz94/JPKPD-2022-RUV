Sat Sep 25 14:41:00 CDT 2021
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
$DATA ../../../../data/spa/D/dat80.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m80.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   21721.2714892518        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.9001E+02  4.0431E+02 -6.9210E+01  3.1810E+02  9.8200E+01 -2.2010E+03 -9.5572E+02 -2.7239E+01 -1.5399E+03 -4.4762E+02
            -4.1345E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -549.565081399014        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.3210E+00  1.1373E+00  8.8889E-01  1.5619E+00  1.4554E+00  1.5197E+00  9.8881E-01  9.4627E-01  8.0100E-01  9.2780E-01
             1.5444E+01
 PARAMETER:  3.7836E-01  2.2862E-01 -1.7779E-02  5.4588E-01  4.7530E-01  5.1849E-01  8.8751E-02  4.4769E-02 -1.2189E-01  2.5057E-02
             2.8372E+00
 GRADIENT:  -2.9540E+01  2.8807E+01 -6.3303E+00  5.5569E+01 -3.8458E+00  1.3634E+01  1.0653E+00  3.9880E+00  4.8560E+00  9.0479E-01
             1.8460E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -560.564602012027        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  1.3232E+00  8.2552E-01  9.9108E-01  1.5411E+00  1.4345E+00  1.3429E+00  9.3099E-01  2.3194E-01  4.7991E-01  8.3531E-01
             1.5972E+01
 PARAMETER:  3.8009E-01 -9.1742E-02  9.1038E-02  5.3250E-01  4.6079E-01  3.9482E-01  2.8490E-02 -1.3613E+00 -6.3415E-01 -7.9953E-02
             2.8708E+00
 GRADIENT:  -4.5435E+00  8.3242E+00 -1.8746E+00  8.5202E+00 -1.0731E+00  5.2509E+00  1.4159E+00  2.1631E-01  3.8809E+00  5.5736E-01
             5.5934E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -568.866786772166        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      227
 NPARAMETR:  1.1866E+00  2.5695E-01  1.2676E+00  1.7377E+00  2.0099E+00  1.2553E+00  1.8368E+00  2.5989E-02  2.2490E-01  3.1358E+00
             1.4332E+01
 PARAMETER:  2.7106E-01 -1.2589E+00  3.3712E-01  6.5257E-01  7.9808E-01  3.2734E-01  7.0804E-01 -3.5501E+00 -1.3921E+00  1.2429E+00
             2.7625E+00
 GRADIENT:  -2.8463E+01  6.1691E+00  3.8797E+00  4.1856E+01 -4.5359E+00 -8.5318E+00  4.2844E-01  7.8593E-04  9.0245E-01  1.6068E+00
            -1.1871E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -630.614233535264        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      303
 NPARAMETR:  7.3541E-01  1.0000E-02  8.7100E-02  7.0415E-01  2.0401E+01  1.4308E+00  4.0834E-02  1.0000E-02  1.0000E-02  4.5131E-01
             1.3152E+01
 PARAMETER: -2.0733E-01 -6.1581E+00 -2.3407E+00 -2.5076E-01  3.1156E+00  4.5825E-01 -3.0982E+00 -2.3237E+01 -8.3858E+00 -6.9559E-01
             2.6766E+00
 GRADIENT:   2.2527E+01  0.0000E+00 -4.7789E+00  2.7759E+01  1.1788E-01  5.9563E+00  7.3267E-07  0.0000E+00  0.0000E+00 -6.4000E-05
            -3.0948E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -641.628963244831        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      375
 NPARAMETR:  4.2089E-01  1.0000E-02  1.9137E-02  2.4500E-01  1.6913E+02  1.3341E+00  1.0000E-02  1.0000E-02  1.0000E-02  2.3499E-01
             1.3241E+01
 PARAMETER: -7.6539E-01 -9.1326E+00 -3.8561E+00 -1.3065E+00  5.2307E+00  3.8823E-01 -7.4923E+00 -3.4717E+01 -1.3665E+01 -1.3482E+00
             2.6833E+00
 GRADIENT:  -6.1385E+00  0.0000E+00 -2.7389E+01  4.8962E+01  2.1997E-03 -1.5475E+00  0.0000E+00  0.0000E+00  0.0000E+00 -4.6799E-08
             7.7142E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -642.302801662148        NO. OF FUNC. EVALS.: 159
 CUMULATIVE NO. OF FUNC. EVALS.:      534
 NPARAMETR:  4.0797E-01  1.0000E-02  1.7592E-02  2.2342E-01  2.1572E+02  1.3298E+00  1.0000E-02  1.0000E-02  1.0000E-02  2.7077E-01
             1.3052E+01
 PARAMETER: -7.9657E-01 -9.2260E+00 -3.9403E+00 -1.3987E+00  5.4740E+00  3.8501E-01 -7.6941E+00 -3.5359E+01 -1.3858E+01 -1.2065E+00
             2.6689E+00
 GRADIENT:   2.7605E-01  0.0000E+00 -5.6675E-01  6.9996E-01  1.1823E-04 -5.8675E-01  0.0000E+00  0.0000E+00  0.0000E+00 -2.1814E-08
            -2.4329E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -642.307898360861        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:      720
 NPARAMETR:  4.0638E-01  1.0000E-02  1.7386E-02  2.2134E-01  2.5642E+02  1.3312E+00  1.0000E-02  1.0000E-02  1.0000E-02  2.6334E-01
             1.3095E+01
 PARAMETER: -8.0046E-01 -9.2437E+00 -3.9521E+00 -1.4081E+00  5.6468E+00  3.8611E-01 -7.7628E+00 -3.5421E+01 -1.3914E+01 -1.2343E+00
             2.6722E+00
 GRADIENT:   1.2673E-01  0.0000E+00 -2.4290E-01  2.2126E-01 -1.2837E-05 -1.5733E-02  0.0000E+00  0.0000E+00  0.0000E+00 -1.2619E-08
            -9.1035E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -642.307936840213        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      895
 NPARAMETR:  4.0633E-01  1.0000E-02  1.7391E-02  2.2135E-01  9.4141E+02  1.3312E+00  1.0000E-02  1.0000E-02  1.0000E-02  2.6334E-01
             1.3095E+01
 PARAMETER: -8.0059E-01 -9.2437E+00 -3.9518E+00 -1.4080E+00  6.9474E+00  3.8612E-01 -7.7628E+00 -3.5421E+01 -1.3914E+01 -1.2343E+00
             2.6722E+00
 GRADIENT:   1.5427E-02  0.0000E+00 -1.0911E-02 -5.4289E-03 -1.2484E-06  2.4551E-04  0.0000E+00  0.0000E+00  0.0000E+00 -9.6712E-10
            -1.0971E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -642.307937172887        NO. OF FUNC. EVALS.: 191
 CUMULATIVE NO. OF FUNC. EVALS.:     1086
 NPARAMETR:  4.0632E-01  1.0000E-02  1.7391E-02  2.2135E-01  1.3127E+03  1.3312E+00  1.0000E-02  1.0000E-02  1.0000E-02  2.6302E-01
             1.3095E+01
 PARAMETER: -8.0060E-01 -9.2437E+00 -3.9518E+00 -1.4080E+00  7.2798E+00  3.8612E-01 -7.7628E+00 -3.5421E+01 -1.3914E+01 -1.2355E+00
             2.6723E+00
 GRADIENT:   1.2940E-02  0.0000E+00 -1.9502E-02  7.1059E-03 -3.7488E-08  7.0129E-04  0.0000E+00  0.0000E+00  0.0000E+00 -5.0608E-10
            -2.6770E-03

0ITERATION NO.:   48    OBJECTIVE VALUE:  -642.307937428226        NO. OF FUNC. EVALS.: 103
 CUMULATIVE NO. OF FUNC. EVALS.:     1189
 NPARAMETR:  4.0631E-01  1.0000E-02  1.7391E-02  2.2134E-01  1.3031E+03  1.3312E+00  1.0000E-02  1.0000E-02  1.0000E-02  2.6287E-01
             1.3095E+01
 PARAMETER: -8.0063E-01 -9.2437E+00 -3.9518E+00 -1.4081E+00  7.2725E+00  3.8611E-01 -7.7628E+00 -3.5421E+01 -1.3914E+01 -1.2361E+00
             2.6723E+00
 GRADIENT:  -5.6881E-03  0.0000E+00  2.5105E-04 -8.3885E-03 -9.3088E-07 -4.9513E-04  0.0000E+00  0.0000E+00  0.0000E+00 -5.0585E-10
             5.9181E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1189
 NO. OF SIG. DIGITS IN FINAL EST.:  4.0
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.9168E-03  3.1256E-06  8.8508E-05 -1.7058E-04  9.5237E-08
 SE:             2.8338E-02  1.5024E-06  2.5921E-04  3.0336E-04  2.1944E-07
 N:                     100         100         100         100         100

 P VAL.:         9.4607E-01  3.7485E-02  7.3276E-01  5.7392E-01  6.6430E-01

 ETASHRINKSD(%)  5.0644E+00  9.9995E+01  9.9132E+01  9.8984E+01  9.9999E+01
 ETASHRINKVR(%)  9.8724E+00  1.0000E+02  9.9992E+01  9.9990E+01  1.0000E+02
 EBVSHRINKSD(%)  5.3705E+00  9.9994E+01  9.9048E+01  9.8866E+01  9.9999E+01
 EBVSHRINKVR(%)  1.0453E+01  1.0000E+02  9.9991E+01  9.9987E+01  1.0000E+02
 RELATIVEINF(%)  3.7172E-01  6.9622E-08  1.6533E-05  1.8430E-05  0.0000E+00
 EPSSHRINKSD(%)  5.7155E+00
 EPSSHRINKVR(%)  1.1104E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -642.30793742822584     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       92.842889135512337     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    14.88
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.43
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -642.308       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         4.06E-01  1.00E-02  1.74E-02  2.21E-01  1.30E+03  1.33E+00  1.00E-02  1.00E-02  1.00E-02  2.63E-01  1.31E+01
 


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
+        3.34E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -2.29E+04  0.00E+00  2.23E+06
 
 TH 4
+       -8.07E+02  0.00E+00 -1.97E+05  2.01E+04
 
 TH 5
+        3.46E-05  0.00E+00 -6.14E-04  3.07E-05 -1.33E-10
 
 TH 6
+       -1.02E+01  0.00E+00  8.75E+02 -1.02E+02  1.24E-06  8.55E+01
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH10
+        1.15E-02  0.00E+00  4.26E-01  3.10E-02 -8.73E-07  4.37E-02  0.00E+00  0.00E+00  0.00E+00 -8.06E-02
 
 TH11
+       -3.03E+01  0.00E+00  4.44E+02 -1.73E+01 -3.66E-07  9.15E-01  0.00E+00  0.00E+00  0.00E+00  7.93E-05  2.22E+00
 
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
 #CPUT: Total CPU Time in Seconds,       21.371
Stop Time:
Sat Sep 25 14:41:31 CDT 2021
