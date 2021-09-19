Sat Sep 18 07:24:45 CDT 2021
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
$DATA ../../../../data/int/D/dat70.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:     1000
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

 TOT. NO. OF OBS RECS:      900
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
 RAW OUTPUT FILE (FILE): m70.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   26494.2579514539        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -8.3812E+01  4.2618E+02 -5.8251E+01  6.6919E+01  2.9197E+02 -3.1245E+03 -1.2912E+03 -6.8194E+01 -2.2564E+03 -1.0700E+03
            -5.3012E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -965.206921092377        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  2.8003E+00  1.6566E+00  9.9178E-01  2.7848E+00  9.7513E-01  5.6697E+00  6.5788E+00  9.6573E-01  6.4936E+00  2.8808E+00
             9.6791E+00
 PARAMETER:  1.1297E+00  6.0479E-01  9.1746E-02  1.1242E+00  7.4811E-02  1.8351E+00  1.9839E+00  6.5129E-02  1.9708E+00  1.1581E+00
             2.3700E+00
 GRADIENT:   5.3622E+01 -2.1183E+00 -4.2200E+01  4.9368E+01 -6.1450E+01  1.1479E+02  8.6488E+01  3.8034E+00  1.1479E+02  6.8556E+01
             2.9687E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1064.51329615472        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      153
 NPARAMETR:  1.4969E+00  6.2679E+00  1.7419E+01  7.5799E+00  5.1579E+00  3.3683E+00  7.3420E+00  5.4542E-01  1.8896E+01  1.9254E+00
             1.0073E+01
 PARAMETER:  5.0341E-01  1.9354E+00  2.9576E+00  2.1255E+00  1.7405E+00  1.3144E+00  2.0936E+00 -5.0620E-01  3.0390E+00  7.5513E-01
             2.4099E+00
 GRADIENT:   1.8525E+01  2.9080E+01 -2.9368E+00  3.4138E+01  5.3323E+01  5.4518E+01  1.0830E+02  1.4924E-02  4.3982E+01  3.5472E+01
             3.0850E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1239.55035081432        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      225
 NPARAMETR:  1.4889E+00  1.5801E+00  5.0407E+00  8.3943E-01  2.3574E+00  3.6834E+00  3.3086E+00  1.6913E-01  4.1201E+00  8.6298E-01
             9.4568E+00
 PARAMETER:  4.9804E-01  5.5746E-01  1.7175E+00 -7.5028E-02  9.5756E-01  1.4038E+00  1.2965E+00 -1.6771E+00  1.5159E+00 -4.7359E-02
             2.3467E+00
 GRADIENT:   1.0022E+01 -8.5548E+01 -1.3147E+01  1.3218E+01  2.4487E+01  5.2372E+01 -1.2747E+02  1.5293E-02  5.6775E+01  1.3464E+01
             3.5379E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1325.23220110852        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      295
 NPARAMETR:  1.0992E+00  1.4762E+00  6.3209E+00  1.0169E+00  2.0431E+00  3.1888E+00  4.7279E+00  1.2019E-01  2.0969E+00  7.6221E-01
             7.5882E+00
 PARAMETER:  1.9456E-01  4.8950E-01  1.9439E+00  1.1671E-01  8.1445E-01  1.2596E+00  1.6535E+00 -2.0187E+00  8.4048E-01 -1.7154E-01
             2.1266E+00
 GRADIENT:  -4.0197E+01 -4.4019E-01  2.1559E+00  4.7283E+00 -1.9885E+01 -1.2327E+01 -1.7895E+01  8.3878E-04  1.5718E+01  1.1557E+01
             1.9711E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1342.31132631356        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      365
 NPARAMETR:  1.2813E+00  1.2623E+00  4.0253E+00  9.6842E-01  1.9628E+00  3.2093E+00  5.3419E+00  1.3557E-01  9.3916E-01  3.0041E-01
             7.7858E+00
 PARAMETER:  3.4787E-01  3.3291E-01  1.4926E+00  6.7912E-02  7.7437E-01  1.2661E+00  1.7756E+00 -1.8983E+00  3.7232E-02 -1.1026E+00
             2.1523E+00
 GRADIENT:  -7.0120E+00 -5.2861E+00 -1.2150E+00 -8.1959E+00 -1.1853E+00 -1.5962E+00  3.9986E+00  5.5909E-03  1.2882E+00  1.7774E+00
             1.4289E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1343.18117909662        NO. OF FUNC. EVALS.:  94
 CUMULATIVE NO. OF FUNC. EVALS.:      459
 NPARAMETR:  1.3183E+00  1.3405E+00  4.3348E+00  9.7306E-01  2.0067E+00  3.2295E+00  5.3052E+00  1.4575E-01  8.6418E-01  2.3530E-01
             7.7727E+00
 PARAMETER:  3.7636E-01  3.9301E-01  1.5667E+00  7.2687E-02  7.9650E-01  1.2723E+00  1.7687E+00 -1.8259E+00 -4.5980E-02 -1.3469E+00
             2.1506E+00
 GRADIENT:  -4.2036E+00 -1.6391E+00  1.3492E-01  7.9345E-01 -2.2050E+00 -9.4422E+00 -1.3053E+01  5.2211E-03 -1.1178E+00  1.0171E+00
            -4.3441E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1344.15402318797        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      637
 NPARAMETR:  1.3477E+00  1.1927E+00  5.6363E+00  1.0537E+00  2.0771E+00  3.2941E+00  5.9004E+00  1.2601E-01  9.4166E-01  1.8075E-01
             7.7890E+00
 PARAMETER:  3.9843E-01  2.7621E-01  1.8292E+00  1.5230E-01  8.3098E-01  1.2921E+00  1.8750E+00 -1.9714E+00  3.9894E-02 -1.6106E+00
             2.1527E+00
 GRADIENT:   2.8508E-01  8.9607E-02  6.0537E-01  1.5380E-01  8.8758E-01 -1.8706E+00 -1.1691E+00  9.2794E-04 -6.0505E-01  5.5634E-01
             2.3308E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1344.46156560954        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      813
 NPARAMETR:  1.3447E+00  1.1892E+00  5.1366E+00  1.0515E+00  2.0366E+00  3.3095E+00  5.9251E+00  4.2769E-02  9.6477E-01  4.4463E-02
             7.7791E+00
 PARAMETER:  3.9620E-01  2.7331E-01  1.7364E+00  1.5023E-01  8.1127E-01  1.2968E+00  1.8792E+00 -3.0519E+00  6.4139E-02 -3.0131E+00
             2.1514E+00
 GRADIENT:  -4.6230E-02 -3.2146E-02 -2.4937E-01 -8.7731E-01  6.2278E-01 -8.0431E-02  7.4537E-01  9.3432E-05  2.4492E-01  3.4279E-02
             1.4627E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1344.48069529454        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      988
 NPARAMETR:  1.3451E+00  1.1838E+00  5.2070E+00  1.0550E+00  2.0395E+00  3.3101E+00  5.9274E+00  1.3149E-02  9.6152E-01  1.0000E-02
             7.7817E+00
 PARAMETER:  3.9650E-01  2.6876E-01  1.7500E+00  1.5357E-01  8.1272E-01  1.2970E+00  1.8796E+00 -4.2314E+00  6.0759E-02 -4.7597E+00
             2.1518E+00
 GRADIENT:  -2.3577E-04 -2.2012E-03 -2.4995E-04  4.8877E-03 -1.2245E-03  2.3472E-03 -5.5276E-03  7.1932E-06 -7.9918E-04  0.0000E+00
            -3.4140E-03

0ITERATION NO.:   46    OBJECTIVE VALUE:  -1344.48069529454        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     1010
 NPARAMETR:  1.3451E+00  1.1838E+00  5.2070E+00  1.0550E+00  2.0395E+00  3.3101E+00  5.9274E+00  1.3149E-02  9.6152E-01  1.0000E-02
             7.7817E+00
 PARAMETER:  3.9650E-01  2.6876E-01  1.7500E+00  1.5357E-01  8.1272E-01  1.2970E+00  1.8796E+00 -4.2314E+00  6.0759E-02 -4.7597E+00
             2.1518E+00
 GRADIENT:  -2.3577E-04 -2.2012E-03 -2.4995E-04  4.8877E-03 -1.2245E-03  2.3472E-03 -5.5276E-03  7.1932E-06 -7.9918E-04  0.0000E+00
            -3.4140E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1010
 NO. OF SIG. DIGITS IN FINAL EST.:  3.5
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         6.0025E-03  2.6209E-02 -5.4121E-05 -6.2922E-02  6.1845E-05
 SE:             2.9443E-02  2.5588E-02  4.1788E-05  1.3297E-02  1.4716E-04
 N:                     100         100         100         100         100

 P VAL.:         8.3846E-01  3.0570E-01  1.9527E-01  2.2253E-06  6.7430E-01

 ETASHRINKSD(%)  1.3608E+00  1.4279E+01  9.9860E+01  5.5453E+01  9.9507E+01
 ETASHRINKVR(%)  2.7030E+00  2.6518E+01  1.0000E+02  8.0156E+01  9.9998E+01
 EBVSHRINKSD(%)  9.3771E-01  8.4134E+00  9.9851E+01  6.2911E+01  9.9446E+01
 EBVSHRINKVR(%)  1.8666E+00  1.6119E+01  1.0000E+02  8.6244E+01  9.9997E+01
 RELATIVEINF(%)  9.8047E+01  4.5093E+01  5.0653E-05  6.9096E+00  7.4992E-04
 EPSSHRINKSD(%)  7.3639E+00
 EPSSHRINKVR(%)  1.4185E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1344.4806952945432     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       309.60866447386752     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    28.17
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    16.67
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1344.481       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.35E+00  1.18E+00  5.21E+00  1.06E+00  2.04E+00  3.31E+00  5.93E+00  1.31E-02  9.62E-01  1.00E-02  7.78E+00
 


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
+        5.44E+01
 
 TH 2
+       -3.95E-01  1.90E+01
 
 TH 3
+        1.04E-01  5.19E-01  7.32E-01
 
 TH 4
+       -3.61E+00  2.60E+01 -2.84E+00  1.99E+02
 
 TH 5
+       -9.59E-01 -6.71E+00 -6.29E+00  1.63E+01  7.58E+01
 
 TH 6
+        6.01E-01 -6.21E-02  1.12E-02 -1.71E-02 -5.65E-01  1.69E+01
 
 TH 7
+        3.47E-02  2.11E+00 -2.34E-01 -1.24E+01  1.96E+00 -7.89E-02  3.60E+00
 
 TH 8
+        2.97E+00 -7.13E+00  5.13E-03 -9.90E+00 -5.80E-01 -5.60E-01 -2.11E-01 -9.78E+00
 
 TH 9
+       -1.87E+00 -2.79E+00 -5.43E-01 -3.48E+01  4.83E+00 -1.59E-01  1.72E+00  2.48E+01  2.25E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -2.65E+00 -2.01E+00 -5.57E-02 -1.29E+01  5.45E-01  8.97E-01  7.12E-01  1.05E-01  5.22E+00  0.00E+00  1.72E+01
 
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
 #CPUT: Total CPU Time in Seconds,       44.959
Stop Time:
Sat Sep 18 07:25:32 CDT 2021
