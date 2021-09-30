Wed Sep 29 20:49:13 CDT 2021
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
$DATA ../../../../data/spa1/B/dat6.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      600
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

 TOT. NO. OF OBS RECS:      500
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
 RAW OUTPUT FILE (FILE): m6.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2131.48153518205        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.3290E+02 -6.3864E+01 -9.9256E+01  6.0412E+01  1.6817E+02  5.1591E+01  1.9833E+00  1.0938E+01  1.3285E+01 -2.4388E+00
             1.2605E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2136.78736043743        NO. OF FUNC. EVALS.: 159
 CUMULATIVE NO. OF FUNC. EVALS.:      172
 NPARAMETR:  1.0597E+00  1.1684E+00  1.2992E+00  9.9549E-01  9.7343E-01  1.0177E+00  9.8055E-01  9.2586E-01  9.7733E-01  9.6553E-01
             9.6643E-01
 PARAMETER:  1.5800E-01  2.5563E-01  3.6178E-01  9.5476E-02  7.3074E-02  1.1757E-01  8.0359E-02  2.2968E-02  7.7071E-02  6.4922E-02
             6.5852E-02
 GRADIENT:   1.2123E+02  9.5990E+01  4.6643E+01  5.4358E+01 -9.5795E+01  1.3980E+00  3.1282E+00 -4.1471E+00 -1.2750E+01 -1.3421E+01
            -2.6966E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2139.85174764959        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      348
 NPARAMETR:  1.0467E+00  1.1153E+00  1.7489E+00  1.0086E+00  1.0855E+00  9.7175E-01  7.0176E-01  9.5927E-01  1.0567E+00  1.2198E+00
             9.9185E-01
 PARAMETER:  1.4569E-01  2.0917E-01  6.5899E-01  1.0854E-01  1.8203E-01  7.1340E-02 -2.5417E-01  5.8413E-02  1.5519E-01  2.9868E-01
             9.1813E-02
 GRADIENT:   1.0413E+02  5.6548E+01  3.9051E+01  3.4050E+01 -5.8151E+01 -1.4468E+01  7.8049E-01 -8.4430E+00 -1.0071E+01  5.4486E+00
            -7.3481E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2147.99195185109        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      524
 NPARAMETR:  9.9966E-01  1.1773E+00  1.0963E+00  9.2721E-01  9.8980E-01  9.9995E-01  7.3388E-01  7.3162E-01  1.1095E+00  1.0487E+00
             9.8499E-01
 PARAMETER:  9.9656E-02  2.6320E-01  1.9194E-01  2.4430E-02  8.9751E-02  9.9954E-02 -2.0941E-01 -2.1250E-01  2.0393E-01  1.4756E-01
             8.4875E-02
 GRADIENT:  -7.1264E+00 -2.4555E+00 -6.7775E+00  3.3444E+00 -5.1634E+00  7.6757E-01 -1.0316E+00  2.2190E+00 -1.4037E+00  2.6740E+00
            -2.4002E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2148.72622416798        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      700
 NPARAMETR:  1.0059E+00  1.4436E+00  9.6850E-01  7.5978E-01  1.0695E+00  9.9833E-01  6.8631E-01  4.7932E-01  1.2947E+00  1.0885E+00
             9.9087E-01
 PARAMETER:  1.0592E-01  4.6716E-01  6.7996E-02 -1.7473E-01  1.6724E-01  9.8325E-02 -2.7643E-01 -6.3539E-01  3.5826E-01  1.8480E-01
             9.0831E-02
 GRADIENT:   3.7637E+00  7.7159E+00 -1.9354E+00  1.0136E+01  2.5843E+00 -4.4701E-01 -4.8315E-01  2.1557E-01 -2.0151E-01  1.3079E+00
             4.2595E-03

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2148.80368133752        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      879             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0057E+00  1.5191E+00  9.0529E-01  6.9851E-01  1.0822E+00  1.0001E+00  6.8920E-01  3.2393E-01  1.3567E+00  1.0793E+00
             9.9144E-01
 PARAMETER:  1.0568E-01  5.1810E-01  4.9867E-04 -2.5881E-01  1.7902E-01  1.0009E-01 -2.7223E-01 -1.0272E+00  4.0502E-01  1.7627E-01
             9.1408E-02
 GRADIENT:   4.8175E+02  4.8268E+02  2.2742E+00  9.2714E+01  1.2104E+01  5.2808E+01  1.0746E+01  9.1855E-02  2.3975E+01  1.2603E+00
             4.5306E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2148.81835477929        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1060             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0061E+00  1.5225E+00  8.9538E-01  6.9751E-01  1.0821E+00  1.0003E+00  6.8831E-01  2.4640E-01  1.3611E+00  1.0809E+00
             9.9205E-01
 PARAMETER:  1.0610E-01  5.2036E-01 -1.0502E-02 -2.6024E-01  1.7895E-01  1.0032E-01 -2.7352E-01 -1.3008E+00  4.0833E-01  1.7780E-01
             9.2020E-02
 GRADIENT:   4.8384E+02  4.8687E+02  1.2241E+00  9.4818E+01  1.4234E+01  5.2872E+01  1.0560E+01  5.2682E-02  2.4586E+01  1.5658E+00
             9.1450E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2148.82169279088        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1243
 NPARAMETR:  1.0054E+00  1.5234E+00  8.9257E-01  6.9821E-01  1.0812E+00  1.0003E+00  6.8858E-01  2.4000E-01  1.3604E+00  1.0811E+00
             9.9215E-01
 PARAMETER:  1.0542E-01  5.2093E-01 -1.3646E-02 -2.5924E-01  1.7811E-01  1.0025E-01 -2.7312E-01 -1.3271E+00  4.0778E-01  1.7796E-01
             9.2116E-02
 GRADIENT:   1.7620E+00 -3.3483E+00  1.3977E-01 -1.5986E-01  7.6950E-01  1.5129E-01 -9.2586E-03  5.6098E-03  2.0285E-01  1.5721E-01
             5.7069E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2148.82370807948        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     1428
 NPARAMETR:  1.0055E+00  1.5228E+00  8.8814E-01  6.9810E-01  1.0798E+00  1.0003E+00  6.9087E-01  2.2007E-01  1.3590E+00  1.0795E+00
             9.9217E-01
 PARAMETER:  1.0553E-01  5.2053E-01 -1.8622E-02 -2.5939E-01  1.7674E-01  1.0031E-01 -2.6980E-01 -1.4138E+00  4.0675E-01  1.7653E-01
             9.2139E-02
 GRADIENT:   1.9850E+00 -4.4556E+00 -2.4244E-01 -4.8519E-01  1.1344E+00  1.7108E-01  4.1758E-02  8.6556E-03  3.2003E-01  2.4222E-01
             1.5928E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2148.82552224590        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1609
 NPARAMETR:  1.0055E+00  1.5230E+00  8.8639E-01  6.9852E-01  1.0784E+00  1.0003E+00  6.9196E-01  1.9411E-01  1.3572E+00  1.0780E+00
             9.9210E-01
 PARAMETER:  1.0553E-01  5.2071E-01 -2.0600E-02 -2.5879E-01  1.7548E-01  1.0031E-01 -2.6823E-01 -1.5393E+00  4.0545E-01  1.7509E-01
             9.2073E-02
 GRADIENT:   1.9696E+00 -2.8316E+00  2.1075E-01  1.7863E-02  1.7268E-01  1.6742E-01 -1.3391E-02  9.6776E-04  1.7012E-01  6.6183E-02
            -1.2765E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -2148.82616502384        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1793             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0055E+00  1.5224E+00  8.8428E-01  6.9862E-01  1.0778E+00  1.0003E+00  6.9302E-01  1.9123E-01  1.3562E+00  1.0773E+00
             9.9212E-01
 PARAMETER:  1.0553E-01  5.2029E-01 -2.2982E-02 -2.5864E-01  1.7490E-01  1.0032E-01 -2.6670E-01 -1.5543E+00  4.0468E-01  1.7442E-01
             9.2087E-02
 GRADIENT:   4.7966E+02  4.8723E+02  5.8187E-01  9.5906E+01  1.3831E+01  5.2908E+01  1.0349E+01  4.2795E-02  2.4265E+01  1.7575E+00
             1.0864E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -2148.82667199074        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     1978             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0056E+00  1.5223E+00  8.8369E-01  6.9865E-01  1.0773E+00  1.0003E+00  6.9361E-01  1.7759E-01  1.3556E+00  1.0767E+00
             9.9210E-01
 PARAMETER:  1.0554E-01  5.2023E-01 -2.3647E-02 -2.5861E-01  1.7450E-01  1.0032E-01 -2.6584E-01 -1.6283E+00  4.0423E-01  1.7392E-01
             9.2070E-02
 GRADIENT:   4.7969E+02  4.8722E+02  7.5885E-01  9.5844E+01  1.3544E+01  5.2911E+01  1.0333E+01  3.6225E-02  2.4190E+01  1.6807E+00
             1.0299E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -2148.82692679870        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     2165
 NPARAMETR:  1.0056E+00  1.5221E+00  8.8299E-01  6.9879E-01  1.0770E+00  1.0003E+00  6.9406E-01  1.7363E-01  1.3550E+00  1.0764E+00
             9.9210E-01
 PARAMETER:  1.0554E-01  5.2010E-01 -2.4437E-02 -2.5841E-01  1.7418E-01  1.0032E-01 -2.6519E-01 -1.6508E+00  4.0383E-01  1.7360E-01
             9.2072E-02
 GRADIENT:   1.9788E+00 -3.5449E+00 -5.4494E-02 -2.1332E-01  4.3130E-01  1.6954E-01  8.1846E-03  2.3682E-03  1.8382E-01  8.6495E-02
             4.2321E-02

0ITERATION NO.:   64    OBJECTIVE VALUE:  -2148.82708697960        NO. OF FUNC. EVALS.: 132
 CUMULATIVE NO. OF FUNC. EVALS.:     2297
 NPARAMETR:  1.0056E+00  1.5220E+00  8.8288E-01  6.9874E-01  1.0767E+00  1.0003E+00  6.9442E-01  1.6324E-01  1.3547E+00  1.0760E+00
             9.9208E-01
 PARAMETER:  1.0554E-01  5.2005E-01 -2.4562E-02 -2.5847E-01  1.7394E-01  1.0033E-01 -2.6468E-01 -1.7126E+00  4.0356E-01  1.7328E-01
             9.2049E-02
 GRADIENT:   1.9870E+00 -3.4495E+00  1.4903E-01 -3.8285E-01  1.4276E-01  1.7129E-01  1.5641E-02  2.6687E-04  1.4137E-01  1.7337E-02
            -1.9186E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     2297
 NO. OF SIG. DIGITS IN FINAL EST.:  2.0

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.1123E-04 -3.3708E-02 -4.0096E-03  1.9810E-02 -2.7909E-02
 SE:             2.9863E-02  1.8619E-02  1.9081E-03  2.4969E-02  2.4866E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9168E-01  7.0228E-02  3.5609E-02  4.2756E-01  2.6171E-01

 ETASHRINKSD(%)  1.0000E-10  3.7624E+01  9.3608E+01  1.6352E+01  1.6695E+01
 ETASHRINKVR(%)  1.0000E-10  6.1093E+01  9.9591E+01  3.0030E+01  3.0603E+01
 EBVSHRINKSD(%)  3.3829E-01  3.6733E+01  9.4598E+01  1.7082E+01  1.4767E+01
 EBVSHRINKVR(%)  6.7543E-01  5.9972E+01  9.9708E+01  3.1247E+01  2.7353E+01
 RELATIVEINF(%)  9.9199E+01  2.5595E+00  5.9520E-02  5.1177E+00  1.8306E+01
 EPSSHRINKSD(%)  3.3202E+01
 EPSSHRINKVR(%)  5.5380E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2148.8270869795956     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1229.8885537749229     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    37.41
 Elapsed covariance  time in seconds:     7.27
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2148.827       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.01E+00  1.52E+00  8.83E-01  6.99E-01  1.08E+00  1.00E+00  6.94E-01  1.63E-01  1.35E+00  1.08E+00  9.92E-01
 


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
 ********************                            STANDARD ERROR OF ESTIMATE                          ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         3.03E-02  5.27E-01  4.12E-01  3.57E-01  1.51E-01  7.09E-02  1.69E-01  2.23E+00  4.82E-01  1.15E-01  4.82E-02
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


         ETA1      ETA2      ETA3      ETA4      ETA5     
 
 ETA1
+       .........
 
 ETA2
+       ......... .........
 
 ETA3
+       ......... ......... .........
 
 ETA4
+       ......... ......... ......... .........
 
 ETA5
+       ......... ......... ......... ......... .........
 


 SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS  ****


         EPS1     
 
 EPS1
+       .........
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


         ETA1      ETA2      ETA3      ETA4      ETA5     
 
 ETA1
+       .........
 
 ETA2
+       ......... .........
 
 ETA3
+       ......... ......... .........
 
 ETA4
+       ......... ......... ......... .........
 
 ETA5
+       ......... ......... ......... ......... .........
 


 SIGMA - CORR MATRIX FOR RANDOM EFFECTS - EPSILONS  ***


         EPS1     
 
 EPS1
+       .........
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                          COVARIANCE MATRIX OF ESTIMATE                         ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        9.19E-04
 
 TH 2
+        1.05E-03  2.78E-01
 
 TH 3
+       -6.83E-04 -1.94E-01  1.70E-01
 
 TH 4
+       -8.18E-04 -1.87E-01  1.32E-01  1.28E-01
 
 TH 5
+        5.61E-05  6.80E-02 -3.59E-02 -4.59E-02  2.27E-02
 
 TH 6
+        1.47E-04 -6.06E-03  4.87E-03  4.37E-03 -1.29E-03  5.03E-03
 
 TH 7
+       -6.79E-04 -2.81E-02  3.82E-03  1.89E-02 -1.10E-02  1.02E-04  2.87E-02
 
 TH 8
+       -5.46E-03 -7.66E-01  7.52E-01  5.25E-01 -1.28E-01  1.94E-02 -5.98E-02  4.96E+00
 
 TH 9
+        1.20E-03  2.40E-01 -1.54E-01 -1.62E-01  6.38E-02 -4.20E-03 -3.85E-02 -5.78E-01  2.32E-01
 
 TH10
+        3.23E-04  3.55E-03  7.34E-03 -2.57E-03  5.80E-03  6.60E-05 -5.12E-03  6.26E-03  7.85E-03  1.33E-02
 
 TH11
+        4.22E-05  9.38E-03 -6.49E-03 -6.32E-03  2.03E-03 -3.61E-04 -8.01E-04 -1.76E-02  7.62E-03 -1.06E-03  2.32E-03
 
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
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                          CORRELATION MATRIX OF ESTIMATE                        ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        3.03E-02
 
 TH 2
+        6.55E-02  5.27E-01
 
 TH 3
+       -5.47E-02 -8.94E-01  4.12E-01
 
 TH 4
+       -7.55E-02 -9.95E-01  8.95E-01  3.57E-01
 
 TH 5
+        1.23E-02  8.57E-01 -5.79E-01 -8.54E-01  1.51E-01
 
 TH 6
+        6.82E-02 -1.62E-01  1.67E-01  1.72E-01 -1.20E-01  7.09E-02
 
 TH 7
+       -1.32E-01 -3.15E-01  5.47E-02  3.13E-01 -4.30E-01  8.46E-03  1.69E-01
 
 TH 8
+       -8.08E-02 -6.53E-01  8.20E-01  6.60E-01 -3.82E-01  1.23E-01 -1.58E-01  2.23E+00
 
 TH 9
+        8.22E-02  9.44E-01 -7.74E-01 -9.41E-01  8.79E-01 -1.23E-01 -4.72E-01 -5.38E-01  4.82E-01
 
 TH10
+        9.22E-02  5.83E-02  1.54E-01 -6.23E-02  3.34E-01  8.06E-03 -2.62E-01  2.44E-02  1.41E-01  1.15E-01
 
 TH11
+        2.89E-02  3.69E-01 -3.27E-01 -3.67E-01  2.79E-01 -1.05E-01 -9.81E-02 -1.64E-01  3.28E-01 -1.90E-01  4.82E-02
 
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
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                      INVERSE COVARIANCE MATRIX OF ESTIMATE                     ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        1.21E+03
 
 TH 2
+        1.81E+01  4.79E+02
 
 TH 3
+       -6.64E+01  1.06E+02  1.79E+02
 
 TH 4
+        1.52E+02  4.93E+02 -1.06E+02  9.39E+02
 
 TH 5
+        1.93E+02 -1.82E+02 -2.38E+02  1.86E+02  6.65E+02
 
 TH 6
+       -4.77E+01 -2.23E+01  2.66E+00 -5.83E+01 -1.10E+01  2.11E+02
 
 TH 7
+        2.38E+01  4.47E+00  2.45E+01 -1.89E+01 -5.20E+01 -9.47E-01  7.56E+01
 
 TH 8
+        3.29E+00 -3.27E+00 -6.79E+00  1.59E-01  2.56E+00  1.06E-01  2.00E+00  8.79E-01
 
 TH 9
+       -2.53E+00 -3.85E+01 -9.73E+00  1.95E+01 -2.49E+01 -9.33E+00  2.99E+01  7.22E-01  6.42E+01
 
 TH10
+       -4.46E+01  1.10E+01 -2.60E+01  7.96E+00 -8.26E+01  4.42E+00  1.37E+01  3.29E+00  4.34E+00  1.30E+02
 
 TH11
+       -2.26E+01 -3.74E+01 -2.80E+01  3.00E+01  3.33E+01  1.52E+01 -6.53E+00 -1.13E+00  8.80E+00  5.30E+01  5.43E+02
 
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
 
 Elapsed finaloutput time in seconds:     0.02
 #CPUT: Total CPU Time in Seconds,       44.753
Stop Time:
Wed Sep 29 20:50:00 CDT 2021
