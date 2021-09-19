Sat Sep 18 10:46:40 CDT 2021
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
$DATA ../../../../data/spa/A3/dat89.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m89.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   389.130789790321        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   8.7826E+01 -1.6146E+01  7.8175E+01 -1.1529E+02  1.5115E+02  4.2207E+01 -2.8537E+01 -3.2079E+01 -9.1704E+01 -1.2490E+02
            -3.7853E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1212.04283776097        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0558E+00  1.0126E+00  9.2673E-01  1.1561E+00  9.7178E-01  7.1730E-01  9.2152E-01  9.8817E-01  9.4862E-01  9.7522E-01
             5.3706E+00
 PARAMETER:  1.5433E-01  1.1253E-01  2.3906E-02  2.4505E-01  7.1376E-02 -2.3227E-01  1.8274E-02  8.8098E-02  4.7249E-02  7.4910E-02
             1.7809E+00
 GRADIENT:   7.4178E+01 -1.8640E+01 -1.8770E+01 -1.0430E+01  6.3153E+00 -1.7048E+01  1.2626E+01  6.4401E+00  2.1364E+01  1.9108E+01
             1.7319E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1235.11033789578        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  1.0391E+00  7.8673E-01  3.6042E-01  1.2108E+00  4.4740E-01  7.9073E-01  3.3099E-01  2.8965E-01  8.9325E-01  3.9951E-01
             5.0096E+00
 PARAMETER:  1.3835E-01 -1.3987E-01 -9.2048E-01  2.9130E-01 -7.0431E-01 -1.3480E-01 -1.0057E+00 -1.1391E+00 -1.2886E-02 -8.1752E-01
             1.7114E+00
 GRADIENT:  -4.7944E+00  6.2798E+01  1.9549E+01  7.4361E+01 -6.6017E+01 -3.2266E+00  7.3725E-01  1.2975E+00  7.7124E+00  7.5345E+00
             1.3265E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1259.24207615820        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      228
 NPARAMETR:  1.0123E+00  6.3050E-01  3.0465E-01  1.1799E+00  3.8899E-01  8.3057E-01  1.9707E-01  1.6993E-01  9.0045E-01  2.3474E-01
             3.9731E+00
 PARAMETER:  1.1219E-01 -3.6125E-01 -1.0886E+00  2.6545E-01 -8.4419E-01 -8.5640E-02 -1.5242E+00 -1.6724E+00 -4.8656E-03 -1.3493E+00
             1.4795E+00
 GRADIENT:   1.4624E+00  1.8377E+01  2.7917E-01  3.1484E+01 -1.9200E+00  3.9635E+00 -2.7683E-01  3.5647E-01 -4.9777E+00  1.2407E+00
            -1.7726E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1263.64911420726        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      299
 NPARAMETR:  1.0098E+00  4.4816E-01  1.7899E-01  1.0710E+00  2.5864E-01  8.2681E-01  4.9925E-02  2.7970E-02  1.0396E+00  1.0535E-01
             4.0401E+00
 PARAMETER:  1.0971E-01 -7.0261E-01 -1.6204E+00  1.6859E-01 -1.2523E+00 -9.0176E-02 -2.8972E+00 -3.4766E+00  1.3880E-01 -2.1504E+00
             1.4963E+00
 GRADIENT:   1.4450E+01 -2.6796E+00  7.0806E+00 -2.5296E+01  7.8382E+00 -5.1621E+00 -1.9508E-02  5.5180E-03  3.5466E+00 -8.8373E-02
             1.3268E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1264.14223491384        NO. OF FUNC. EVALS.: 137
 CUMULATIVE NO. OF FUNC. EVALS.:      436
 NPARAMETR:  1.0066E+00  4.2289E-01  1.7400E-01  1.0902E+00  2.4945E-01  8.4618E-01  4.8155E-02  2.5439E-02  1.0356E+00  1.0192E-01
             3.9829E+00
 PARAMETER:  1.0662E-01 -7.6064E-01 -1.6487E+00  1.8639E-01 -1.2885E+00 -6.7020E-02 -2.9333E+00 -3.5715E+00  1.3498E-01 -2.1836E+00
             1.4820E+00
 GRADIENT:   4.2397E-01  3.5209E-01  2.5298E-01  4.6675E-01 -2.7692E-01 -9.2033E-02 -2.1818E-02  3.2813E-03 -1.4426E-01 -2.2667E-01
            -1.8222E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1264.36469676476        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      616
 NPARAMETR:  1.0006E+00  4.4539E-01  1.5476E-01  1.0537E+00  2.4133E-01  8.5364E-01  3.1495E-02  1.7517E-02  1.0906E+00  2.4498E-01
             3.9152E+00
 PARAMETER:  1.0060E-01 -7.0881E-01 -1.7659E+00  1.5229E-01 -1.3216E+00 -5.8245E-02 -3.3579E+00 -3.9446E+00  1.8673E-01 -1.3066E+00
             1.4649E+00
 GRADIENT:  -3.8728E+00 -2.5312E+00 -3.8773E+00  3.5298E-01  8.3894E+00  3.9787E-01 -4.3450E-03  1.8687E-03  3.8703E-01  1.7323E-01
             2.8662E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1264.47857345141        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      794            RESET HESSIAN, TYPE II
 NPARAMETR:  1.0044E+00  4.2322E-01  1.5707E-01  1.0614E+00  2.3658E-01  8.5354E-01  1.8007E-01  1.0656E-02  1.0854E+00  2.4912E-01
             3.9009E+00
 PARAMETER:  1.0440E-01 -7.5986E-01 -1.7511E+00  1.5957E-01 -1.3415E+00 -5.8365E-02 -1.6144E+00 -4.4416E+00  1.8190E-01 -1.2898E+00
             1.4612E+00
 GRADIENT:   1.6669E+00 -1.0869E+00 -3.0051E-02  3.2428E-01  7.3453E+00  1.0838E-01 -7.7939E-02  7.6844E-04  2.8891E-01  3.2890E-01
             2.8296E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1264.59865141887        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      865
 NPARAMETR:  1.0035E+00  4.2283E-01  1.5605E-01  1.0591E+00  2.3544E-01  8.5414E-01  4.2404E-01  1.0000E-02  1.0855E+00  1.9879E-01
             3.8707E+00
 PARAMETER:  1.0348E-01 -7.6078E-01 -1.7576E+00  1.5742E-01 -1.3463E+00 -5.7658E-02 -7.5792E-01 -4.5920E+00  1.8200E-01 -1.5155E+00
             1.4534E+00
 GRADIENT:   3.0842E-01 -3.7626E-02  8.3712E-01  3.4664E-01  3.3430E+00 -6.8457E-03  1.4921E-02  0.0000E+00 -1.6075E-01  1.0458E-01
            -1.9625E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1264.61288557516        NO. OF FUNC. EVALS.: 137
 CUMULATIVE NO. OF FUNC. EVALS.:     1002
 NPARAMETR:  1.0043E+00  4.2351E-01  1.5615E-01  1.0599E+00  2.3580E-01  8.5388E-01  4.1590E-01  1.0000E-02  1.0837E+00  1.9084E-01
             3.8960E+00
 PARAMETER:  1.0429E-01 -7.5917E-01 -1.7569E+00  1.5814E-01 -1.3448E+00 -5.7961E-02 -7.7732E-01 -4.6079E+00  1.8040E-01 -1.5563E+00
             1.4599E+00
 GRADIENT:  -7.0473E-02 -5.9274E-01 -8.2804E-01  3.4521E-01 -6.0878E-01 -2.3965E-02  1.2823E-02  0.0000E+00 -5.2666E-02  1.0516E-01
             9.0711E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1264.67657563966        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1180
 NPARAMETR:  1.0046E+00  4.4158E-01  1.6521E-01  1.0713E+00  2.4643E-01  8.4969E-01  4.7262E-01  1.0000E-02  1.0564E+00  1.2531E-01
             3.9132E+00
 PARAMETER:  1.0458E-01 -7.1739E-01 -1.7006E+00  1.6887E-01 -1.3007E+00 -6.2882E-02 -6.4947E-01 -4.5833E+00  1.5484E-01 -1.9770E+00
             1.4644E+00
 GRADIENT:   5.3117E-02 -3.6482E-01 -1.3811E-01 -4.5895E-02  6.8302E-01  3.9065E-02  1.6964E-02  0.0000E+00  7.2651E-02  2.7242E-02
             1.9799E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1264.67872118120        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1355
 NPARAMETR:  1.0044E+00  4.4465E-01  1.6539E-01  1.0709E+00  2.4737E-01  8.4924E-01  4.8500E-01  1.0000E-02  1.0547E+00  1.0101E-01
             3.9163E+00
 PARAMETER:  1.0436E-01 -7.1046E-01 -1.6995E+00  1.6848E-01 -1.2969E+00 -6.3419E-02 -6.2360E-01 -4.7831E+00  1.5330E-01 -2.1925E+00
             1.4651E+00
 GRADIENT:   6.3655E-03 -1.4486E-02  1.7232E-02  9.1302E-03 -1.5443E-02  6.0586E-05 -2.1982E-04  0.0000E+00 -1.1292E-03 -7.6536E-04
            -3.2388E-03

0ITERATION NO.:   58    OBJECTIVE VALUE:  -1264.67872315428        NO. OF FUNC. EVALS.:  91
 CUMULATIVE NO. OF FUNC. EVALS.:     1446
 NPARAMETR:  1.0044E+00  4.4466E-01  1.6538E-01  1.0709E+00  2.4737E-01  8.4924E-01  4.8502E-01  1.0000E-02  1.0547E+00  1.0136E-01
             3.9163E+00
 PARAMETER:  1.0436E-01 -7.1044E-01 -1.6995E+00  1.6847E-01 -1.2969E+00 -6.3417E-02 -6.2356E-01 -4.7831E+00  1.5330E-01 -2.1891E+00
             1.4652E+00
 GRADIENT:   2.9552E-03 -4.1538E-02 -9.7187E-03  5.2929E-03  5.2295E-02  5.5621E-04  2.8563E-03  0.0000E+00  5.7742E-04  8.2410E-06
             2.2867E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1446
 NO. OF SIG. DIGITS IN FINAL EST.:  3.0
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         5.3841E-04 -1.7289E-03  1.1093E-04 -1.9740E-02  2.5366E-03
 SE:             2.8166E-02  8.7914E-03  1.7412E-04  2.3783E-02  3.9928E-03
 N:                     100         100         100         100         100

 P VAL.:         9.8475E-01  8.4410E-01  5.2405E-01  4.0655E-01  5.2523E-01

 ETASHRINKSD(%)  5.6409E+00  7.0548E+01  9.9417E+01  2.0323E+01  8.6624E+01
 ETASHRINKVR(%)  1.0964E+01  9.1326E+01  9.9997E+01  3.6515E+01  9.8211E+01
 EBVSHRINKSD(%)  5.5500E+00  7.0730E+01  9.9394E+01  1.9408E+01  8.6724E+01
 EBVSHRINKVR(%)  1.0792E+01  9.1433E+01  9.9996E+01  3.5049E+01  9.8237E+01
 RELATIVEINF(%)  7.8339E+01  2.5709E-01  1.9311E-04  2.0723E+01  2.7819E-02
 EPSSHRINKSD(%)  2.0788E+01
 EPSSHRINKVR(%)  3.7255E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1264.6787231542767     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -529.52789659053849     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    17.56
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.65
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1264.679       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.00E+00  4.45E-01  1.65E-01  1.07E+00  2.47E-01  8.49E-01  4.85E-01  1.00E-02  1.05E+00  1.01E-01  3.92E+00
 


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
+        1.37E+03
 
 TH 2
+       -1.52E+02  1.89E+03
 
 TH 3
+       -6.81E+02  4.44E+03  1.44E+04
 
 TH 4
+       -7.14E+01  1.89E+02 -7.40E+02  5.51E+02
 
 TH 5
+        9.03E+02 -7.45E+03 -1.91E+04 -3.12E+02  3.08E+04
 
 TH 6
+       -5.21E+00 -7.34E+00  5.72E+01 -1.70E+01  6.14E+01  2.22E+02
 
 TH 7
+       -3.33E+00 -5.49E+01 -1.47E+02 -3.70E+00  2.36E+02  5.45E-01  1.07E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        8.71E+00 -3.15E+01  1.17E+02 -2.05E+01  1.39E+02 -1.57E+00  3.06E+00  0.00E+00  7.27E+01
 
 TH10
+       -3.54E+00 -5.68E+01 -1.48E+02 -4.31E+00  2.62E+02 -1.40E+00  6.52E+00  0.00E+00  5.36E-01  9.44E+00
 
 TH11
+       -2.13E+01 -1.68E+01 -4.34E+01 -6.81E+00  8.03E+01  4.28E+00  6.68E+00  0.00E+00  8.53E+00  7.45E+00  2.77E+01
 
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
 #CPUT: Total CPU Time in Seconds,       24.277
Stop Time:
Sat Sep 18 10:47:06 CDT 2021
