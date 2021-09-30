Thu Sep 30 09:12:13 CDT 2021
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
$DATA ../../../../data/spa2/D/dat50.csv ignore=@
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
Current Date:       30 SEP 2021
Days until program expires : 199
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
 NO. OF DATA RECS IN DATA SET:      700
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

 TOT. NO. OF OBS RECS:      600
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
 RAW OUTPUT FILE (FILE): m50.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   17228.3401450554        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.3341E+02  2.7357E+02 -1.2655E+01  1.4817E+02  2.2017E+02 -1.4070E+03 -6.3911E+02 -3.9537E+01 -1.1129E+03 -6.4014E+02
            -3.5248E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -692.279350119446        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.6153E+00  1.2278E+00  9.0221E-01  1.6480E+00  1.0836E+00  2.4215E+00  1.7589E+00  9.7939E-01  1.5464E+00  1.3256E+00
             1.3741E+01
 PARAMETER:  5.7953E-01  3.0520E-01 -2.9080E-03  5.9955E-01  1.8033E-01  9.8438E-01  6.6471E-01  7.9179E-02  5.3593E-01  3.8183E-01
             2.7204E+00
 GRADIENT:   5.9483E+01 -3.2264E+01 -3.3078E+01  2.1727E+01  4.6391E+01  7.7919E+01 -1.2973E+01  4.3808E+00 -1.7515E+01  2.1075E+01
             3.0748E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -764.920915214007        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      159
 NPARAMETR:  1.5517E+00  1.6966E+00  1.0411E+01  1.6691E+00  4.4080E+00  2.5494E+00  4.1543E+00  5.1919E-01  4.7188E+00  3.5874E+00
             1.1778E+01
 PARAMETER:  5.3932E-01  6.2865E-01  2.4429E+00  6.1227E-01  1.5834E+00  1.0359E+00  1.5241E+00 -5.5548E-01  1.6515E+00  1.3774E+00
             2.5663E+00
 GRADIENT:   5.5246E+01  2.4232E+00 -2.1469E+00  1.2110E+01 -3.7898E+01  5.3688E+01  4.6390E+01  6.3424E-02  7.6676E+01  5.2059E+01
             2.8181E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -861.051153964032        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      231
 NPARAMETR:  1.0362E+00  2.5206E+00  2.1391E+01  7.8931E-01  6.7384E+00  2.1411E+00  1.7726E+00  5.1603E-01  4.3942E+00  3.0718E+00
             9.5823E+00
 PARAMETER:  1.3555E-01  1.0245E+00  3.1630E+00 -1.3660E-01  2.0078E+00  8.6130E-01  6.7246E-01 -5.6160E-01  1.5803E+00  1.2223E+00
             2.3599E+00
 GRADIENT:  -5.9832E+01  5.6391E+01 -1.8416E+00  7.6131E-01 -2.3450E+00  4.8679E+01  1.5883E+01  1.5962E-02 -1.8756E+01  3.2560E+00
             1.5989E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -878.362149573432        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      302
 NPARAMETR:  1.0984E+00  1.7018E+00  4.5099E+01  1.0274E+00  8.5987E+00  1.6797E+00  1.3943E+00  5.2240E-01  4.0046E+00  2.9714E+00
             8.5192E+00
 PARAMETER:  1.9381E-01  6.3171E-01  3.9089E+00  1.2705E-01  2.2516E+00  6.1861E-01  4.3237E-01 -5.4932E-01  1.4874E+00  1.1890E+00
             2.2423E+00
 GRADIENT:  -1.0714E+01 -6.9842E-01  2.5812E-01 -8.5672E+00  1.4242E+00 -5.3837E+00  1.6910E+01  1.6782E-04 -2.4321E+00  1.8810E-01
             2.1524E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -879.214984781934        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      376
 NPARAMETR:  1.1087E+00  1.5772E+00  5.7558E+01  1.1218E+00  9.0573E+00  1.6910E+00  1.1605E+00  1.2139E+00  3.8725E+00  2.7654E+00
             8.3597E+00
 PARAMETER:  2.0317E-01  5.5563E-01  4.1528E+00  2.1495E-01  2.3036E+00  6.2533E-01  2.4886E-01  2.9385E-01  1.4539E+00  1.1172E+00
             2.2234E+00
 GRADIENT:  -4.1905E-01 -1.2401E+00  2.7398E-01 -5.0145E+00  1.4714E+00 -3.7688E+00  1.2143E+01 -8.6271E-05 -4.9991E+00  6.4078E-02
            -5.9444E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -880.284164229244        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      451
 NPARAMETR:  1.1144E+00  1.4324E+00  9.1552E+01  1.2451E+00  9.8853E+00  1.7120E+00  7.6783E-01  1.1869E+01  3.7592E+00  2.3920E+00
             8.2101E+00
 PARAMETER:  2.0834E-01  4.5935E-01  4.6169E+00  3.1922E-01  2.3910E+00  6.3766E-01 -1.6419E-01  2.5739E+00  1.4242E+00  9.7212E-01
             2.2054E+00
 GRADIENT:   6.7676E+00 -1.6308E+00  1.2917E-01  1.2485E+00  1.1438E+00 -4.1552E-01  5.9583E+00  1.9962E-01 -6.4430E+00  5.2126E-03
            -3.2737E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -881.607946669318        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      525
 NPARAMETR:  1.1202E+00  1.4404E+00  1.1247E+02  1.2587E+00  1.0792E+01  1.6905E+00  6.2813E-01  1.3073E+01  3.7666E+00  2.3329E+00
             8.2825E+00
 PARAMETER:  2.1348E-01  4.6489E-01  4.8227E+00  3.3009E-01  2.4788E+00  6.2503E-01 -3.6501E-01  2.6706E+00  1.4262E+00  9.4712E-01
             2.2141E+00
 GRADIENT:   8.8435E+00  3.4996E+00 -6.4613E-02  3.8570E+00  7.9098E-01 -4.2937E+00  3.9293E+00  5.5650E-01 -9.8749E+00 -6.8167E-03
            -2.4609E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -882.119572265510        NO. OF FUNC. EVALS.: 126
 CUMULATIVE NO. OF FUNC. EVALS.:      651
 NPARAMETR:  1.1191E+00  1.4407E+00  1.1514E+02  1.2549E+00  1.0886E+01  1.6912E+00  5.6880E-01  1.3212E+01  3.7735E+00  3.7610E+00
             8.2984E+00
 PARAMETER:  2.1255E-01  4.6513E-01  4.8461E+00  3.2707E-01  2.4875E+00  6.2545E-01 -4.6423E-01  2.6811E+00  1.4280E+00  1.4247E+00
             2.2161E+00
 GRADIENT:   7.7552E+00  4.4299E+00 -4.0850E-01  3.7327E+00  7.6119E-01 -4.1786E+00  3.1807E+00  1.2988E+00 -1.0893E+01 -1.1707E-02
            -2.2168E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -884.097483936748        NO. OF FUNC. EVALS.:  78
 CUMULATIVE NO. OF FUNC. EVALS.:      729
 NPARAMETR:  1.1079E+00  1.3933E+00  1.1473E+02  1.2530E+00  1.0952E+01  1.6936E+00  1.4839E-01  1.3108E+01  3.8024E+00  2.6876E+00
             8.4176E+00
 PARAMETER:  2.0244E-01  4.3171E-01  4.8426E+00  3.2553E-01  2.4935E+00  6.2685E-01 -1.8079E+00  2.6733E+00  1.4356E+00  1.0887E+00
             2.2303E+00
 GRADIENT:  -1.0107E+00  7.6773E-01  6.9373E-02  7.0653E-01  5.8728E-01 -2.6573E+00  2.0548E-01  2.7156E-01 -8.0501E+00 -1.1911E-02
             1.2006E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -884.848739978075        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      801
 NPARAMETR:  1.1140E+00  1.3944E+00  1.1425E+02  1.2402E+00  1.0994E+01  1.7008E+00  4.0297E-02  1.2991E+01  3.8505E+00  5.2185E+00
             8.4439E+00
 PARAMETER:  2.0792E-01  4.3246E-01  4.8384E+00  3.1531E-01  2.4974E+00  6.3108E-01 -3.1115E+00  2.6643E+00  1.4482E+00  1.7522E+00
             2.2334E+00
 GRADIENT:   2.6372E+00 -2.1914E+00  1.3766E-01 -4.1769E-01  3.5105E-01 -1.0898E+00  1.7761E-02  1.1560E-01 -3.8196E+00  6.3928E-02
             5.9330E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -890.524686871902        NO. OF FUNC. EVALS.: 141
 CUMULATIVE NO. OF FUNC. EVALS.:      942
 NPARAMETR:  1.1474E+00  1.5962E+00  1.0715E+02  1.0634E+00  1.0993E+01  1.8273E+00  1.0000E-02  1.1796E+01  4.7467E+00  1.5907E+00
             8.6370E+00
 PARAMETER:  2.3751E-01  5.6764E-01  4.7743E+00  1.6146E-01  2.4972E+00  7.0284E-01 -2.2148E+01  2.5677E+00  1.6574E+00  5.6416E-01
             2.2561E+00
 GRADIENT:   4.4850E+00 -4.7000E+00  1.2878E-01  5.3574E-01  6.8367E-01  1.1707E+01  0.0000E+00 -1.7746E-03 -5.2648E+00 -4.8216E-03
             3.7515E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -891.228245581276        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1119
 NPARAMETR:  1.1427E+00  1.6707E+00  1.0912E+02  9.5620E-01  1.0568E+01  1.7423E+00  1.0000E-02  1.1314E+01  5.0787E+00  1.5672E+00
             8.6869E+00
 PARAMETER:  2.3336E-01  6.1325E-01  4.7925E+00  5.5212E-02  2.4578E+00  6.5522E-01 -2.7086E+01  2.5260E+00  1.7251E+00  5.4930E-01
             2.2618E+00
 GRADIENT:   2.8288E+00 -2.1012E+00  1.2934E-01 -1.1114E+00  8.0702E-01 -1.1892E+00  0.0000E+00 -2.5396E-03  8.3026E-01 -6.0046E-03
             3.6682E+00

0ITERATION NO.:   65    OBJECTIVE VALUE:  -891.238992132310        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1294
 NPARAMETR:  1.1368E+00  1.6642E+00  1.0913E+02  9.7480E-01  1.0564E+01  1.7477E+00  1.0000E-02  1.1357E+01  5.0465E+00  1.5967E+00
             8.6593E+00
 PARAMETER:  2.2823E-01  6.0937E-01  4.7926E+00  7.4473E-02  2.4575E+00  6.5830E-01 -2.6538E+01  2.5298E+00  1.7187E+00  5.6793E-01
             2.2586E+00
 GRADIENT:  -6.2793E-03 -1.2063E-01  1.2893E-01 -7.0192E-02  8.1343E-01 -3.6886E-01  0.0000E+00 -2.5023E-03  1.9860E-01 -6.1515E-03
            -2.6597E-01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -891.243513684672        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1473
 NPARAMETR:  1.1339E+00  1.6619E+00  1.0780E+02  9.8775E-01  1.0396E+01  1.7474E+00  1.0000E-02  1.1368E+01  5.0376E+00  2.0736E+00
             8.6431E+00
 PARAMETER:  2.2562E-01  6.0795E-01  4.7802E+00  8.7674E-02  2.4414E+00  6.5813E-01 -2.6241E+01  2.5308E+00  1.7169E+00  8.2927E-01
             2.2568E+00
 GRADIENT:  -1.4535E+00  1.3041E+00  1.2840E-01  8.2385E-01  8.8158E-01 -4.7372E-01  0.0000E+00 -1.8600E-03  2.4008E-01 -9.1230E-03
            -2.9532E+00

0ITERATION NO.:   75    OBJECTIVE VALUE:  -891.320298403524        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1651
 NPARAMETR:  1.1375E+00  1.6673E+00  1.0196E+02  9.7181E-01  9.7905E+00  1.7493E+00  1.0000E-02  1.1341E+01  5.0574E+00  4.1761E+00
             8.6585E+00
 PARAMETER:  2.2881E-01  6.1118E-01  4.7246E+00  7.1406E-02  2.3814E+00  6.5922E-01 -2.6198E+01  2.5285E+00  1.7208E+00  1.5294E+00
             2.2585E+00
 GRADIENT:   1.8120E-02 -1.4795E-01  1.2323E-01 -1.2516E-01  7.2235E-01 -2.1872E-02  0.0000E+00  2.7283E-03  4.8664E-01  8.8980E-02
            -3.2003E-01

0ITERATION NO.:   80    OBJECTIVE VALUE:  -900.747113989783        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1832
 NPARAMETR:  1.1065E+00  1.6103E+00  3.5256E+01  1.0029E+00  2.5995E+00  1.7135E+00  1.0000E-02  1.1630E+01  5.1995E+00  1.5272E-01
             8.3155E+00
 PARAMETER:  2.0118E-01  5.7640E-01  3.6626E+00  1.0289E-01  1.0553E+00  6.3855E-01 -1.9877E+01  2.5536E+00  1.7486E+00 -1.7791E+00
             2.2181E+00
 GRADIENT:  -1.5128E+01  1.9762E+01 -4.9160E+00  7.5272E+00  1.7519E+00 -8.5806E+00  0.0000E+00  8.8875E+00  8.6464E+00  1.0781E-01
            -3.0769E+01

0ITERATION NO.:   85    OBJECTIVE VALUE:  -902.270958223210        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     2013
 NPARAMETR:  1.1280E+00  1.5929E+00  3.0509E+01  9.6105E-01  2.4809E+00  1.7029E+00  1.0000E-02  9.7205E+00  4.7957E+00  1.0000E-02
             8.4997E+00
 PARAMETER:  2.2045E-01  5.6557E-01  3.5180E+00  6.0268E-02  1.0086E+00  6.3230E-01 -2.1204E+01  2.3742E+00  1.6677E+00 -2.2281E+01
             2.2400E+00
 GRADIENT:  -6.1862E+00  2.2006E+01 -5.0365E+00  8.1863E-01 -1.7457E+00 -8.6656E+00  0.0000E+00  8.1686E+00 -1.0987E+01  0.0000E+00
            -2.2096E+00

0ITERATION NO.:   90    OBJECTIVE VALUE:  -909.402368010544        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     2193
 NPARAMETR:  1.1217E+00  1.5507E+00  1.6607E+01  9.4819E-01  2.3364E+00  1.7171E+00  1.0000E-02  1.1726E+00  5.0593E+00  1.0000E-02
             8.3998E+00
 PARAMETER:  2.1484E-01  5.3869E-01  2.9098E+00  4.6800E-02  9.4862E-01  6.4064E-01 -1.0579E+02  2.5923E-01  1.7212E+00 -2.1679E+02
             2.2282E+00
 GRADIENT:  -3.7151E+00 -1.0524E+00 -7.1342E+00 -9.4087E-01  7.5577E+00 -6.5940E+00  0.0000E+00  1.5970E+00  8.7692E-01  0.0000E+00
            -3.1503E+00

0ITERATION NO.:   95    OBJECTIVE VALUE:  -912.290988313877        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     2372
 NPARAMETR:  1.1296E+00  1.5294E+00  3.1792E+01  9.7478E-01  2.2671E+00  1.7579E+00  1.0000E-02  1.3878E-01  4.9718E+00  1.0000E-02
             8.4011E+00
 PARAMETER:  2.2187E-01  5.2484E-01  3.5592E+00  7.4457E-02  9.1849E-01  6.6409E-01 -2.1783E+02 -1.8748E+00  1.7038E+00 -5.1633E+02
             2.2284E+00
 GRADIENT:   4.0882E-02 -2.6421E-02 -1.4421E-01  1.6079E-02  8.9238E-02  5.0659E-02  0.0000E+00  1.2151E-02  1.1057E-01  0.0000E+00
            -3.7504E-01

0ITERATION NO.:  100    OBJECTIVE VALUE:  -912.305714742918        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     2554
 NPARAMETR:  1.1305E+00  1.5336E+00  3.2661E+01  9.7357E-01  2.2661E+00  1.7604E+00  1.0000E-02  4.1035E-02  4.9993E+00  1.0000E-02
             8.4060E+00
 PARAMETER:  2.2270E-01  5.2762E-01  3.5862E+00  7.3215E-02  9.1806E-01  6.6556E-01 -2.2131E+02 -3.0933E+00  1.7093E+00 -5.2587E+02
             2.2289E+00
 GRADIENT:   4.1849E-01  5.2464E-01  9.5723E-03  3.4118E-01 -3.3381E-01  5.1251E-01  0.0000E+00  1.0158E-03  1.0219E+00  0.0000E+00
            -1.3069E-01

0ITERATION NO.:  105    OBJECTIVE VALUE:  -912.310273287023        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     2735
 NPARAMETR:  1.1305E+00  1.5361E+00  3.2645E+01  9.6833E-01  2.2683E+00  1.7597E+00  1.0000E-02  1.7668E-02  5.0180E+00  1.0000E-02
             8.4081E+00
 PARAMETER:  2.2267E-01  5.2922E-01  3.5857E+00  6.7815E-02  9.1902E-01  6.6513E-01 -2.2131E+02 -3.9360E+00  1.7130E+00 -5.2587E+02
             2.2292E+00
 GRADIENT:   4.2038E-01  3.1430E-02  8.4247E-03  1.5909E-01 -1.0656E-01  3.9662E-01  0.0000E+00  1.8865E-04  1.5946E+00  0.0000E+00
             2.9715E-01

0ITERATION NO.:  110    OBJECTIVE VALUE:  -912.312361396297        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     2920
 NPARAMETR:  1.1295E+00  1.5376E+00  3.2550E+01  9.6369E-01  2.2685E+00  1.7579E+00  1.0000E-02  1.0000E-02  5.0377E+00  1.0000E-02
             8.4021E+00
 PARAMETER:  2.2174E-01  5.3025E-01  3.5828E+00  6.3010E-02  9.1911E-01  6.6411E-01 -2.2131E+02 -4.5843E+00  1.7169E+00 -5.2587E+02
             2.2285E+00
 GRADIENT:   6.0404E-02 -3.9358E-01  3.2917E-03  5.0752E-02  1.7193E-02  7.0321E-02  0.0000E+00  0.0000E+00  2.2445E+00  0.0000E+00
            -4.2527E-01

0ITERATION NO.:  111    OBJECTIVE VALUE:  -912.312361396297        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     2942
 NPARAMETR:  1.1295E+00  1.5376E+00  3.2550E+01  9.6369E-01  2.2685E+00  1.7579E+00  1.0000E-02  1.0000E-02  5.0377E+00  1.0000E-02
             8.4021E+00
 PARAMETER:  2.2174E-01  5.3025E-01  3.5828E+00  6.3010E-02  9.1911E-01  6.6411E-01 -2.2131E+02 -4.5843E+00  1.7169E+00 -5.2587E+02
             2.2285E+00
 GRADIENT:   6.0404E-02 -3.9358E-01  3.2917E-03  5.0752E-02  1.7193E-02  7.0321E-02  0.0000E+00  0.0000E+00  2.2445E+00  0.0000E+00
            -4.2527E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     2942
 NO. OF SIG. DIGITS IN FINAL EST.:  2.4
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.5053E-02 -9.4252E-04  7.6444E-06  1.3674E-02  5.6240E-05
 SE:             2.7838E-02  2.2552E-04  7.3942E-06  2.6398E-02  7.2823E-05
 N:                     100         100         100         100         100

 P VAL.:         5.8869E-01  2.9250E-05  3.0121E-01  6.0446E-01  4.3994E-01

 ETASHRINKSD(%)  6.7379E+00  9.9244E+01  9.9975E+01  1.1565E+01  9.9756E+01
 ETASHRINKVR(%)  1.3022E+01  9.9994E+01  1.0000E+02  2.1792E+01  9.9999E+01
 EBVSHRINKSD(%)  6.9611E+00  9.9556E+01  9.9940E+01  7.5157E+00  9.9720E+01
 EBVSHRINKVR(%)  1.3438E+01  9.9998E+01  1.0000E+02  1.4467E+01  9.9999E+01
 RELATIVEINF(%)  8.5708E+01  1.0644E-03  3.4539E-05  4.8230E+01  6.8796E-04
 EPSSHRINKSD(%)  8.2306E+00
 EPSSHRINKVR(%)  1.5784E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          600
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1102.7262398456071     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -912.31236139629743     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       190.41387844930966     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    62.64
0R MATRIX ALGORITHMICALLY SINGULAR
0COVARIANCE MATRIX UNOBTAINABLE
0R MATRIX IS OUTPUT
0S MATRIX ALGORITHMICALLY SINGULAR
0S MATRIX IS OUTPUT
0T MATRIX - EQUAL TO RS*RMAT, WHERE S* IS A PSEUDO INVERSE OF S - IS OUTPUT
 Elapsed covariance  time in seconds:    11.68
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -912.312       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.13E+00  1.54E+00  3.25E+01  9.64E-01  2.27E+00  1.76E+00  1.00E-02  1.00E-02  5.04E+00  1.00E-02  8.40E+00
 


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
 ********************                                     T MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        2.93E+02
 
 TH 2
+       -1.00E+02  1.52E+02
 
 TH 3
+       -9.25E-03  8.78E-03  6.34E-07
 
 TH 4
+       -2.57E+01  3.95E+01  2.07E-03  1.09E+01
 
 TH 5
+        8.30E+00 -1.91E+01 -1.06E-03 -4.92E+00  2.49E+00
 
 TH 6
+       -6.32E+01  3.44E+01  8.35E-04  1.38E+01 -3.10E+00  5.71E+01
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        9.32E-30 -4.97E-30 -1.41E-34 -1.95E-30  4.50E-31 -7.89E-30  0.00E+00 -7.34E-56
 
 TH 9
+        9.64E+00 -1.36E+01 -8.46E-04 -3.38E+00  1.71E+00 -1.96E+00  0.00E+00 -1.99E-26  1.25E+00
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -3.09E+00 -8.93E+00 -4.35E-04 -2.11E+00  1.27E+00  9.28E-01  0.00E+00  8.50E-27  8.10E-01  0.00E+00  1.34E+00
 
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
 ********************                                     R MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        2.35E+02
 
 TH 2
+       -1.03E+01  1.66E+02
 
 TH 3
+       -8.92E-03  1.68E-03  2.18E-03
 
 TH 4
+       -4.46E+00  3.76E+01  7.76E-03  4.06E+01
 
 TH 5
+       -1.91E+00 -2.07E+01 -4.83E-02 -5.37E+00  1.79E+01
 
 TH 6
+       -4.29E+00  5.13E+00 -1.68E-03 -3.31E-01 -6.08E-01  4.88E+01
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  1.75E-01
 
 TH 9
+        1.45E+00 -1.55E+01  7.78E-03  2.46E+00  1.56E+00 -6.24E-01  0.00E+00  0.00E+00  5.29E+00
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -7.60E+00 -1.41E+01 -7.17E-04 -3.67E+00  1.24E+00  1.18E+00  0.00E+00  0.00E+00  8.81E-01  0.00E+00  1.04E+01
 
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
 ********************                                     S MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        2.45E+02
 
 TH 2
+        6.86E+01  2.30E+02
 
 TH 3
+       -3.84E-03 -6.47E-04  6.19E-04
 
 TH 4
+        6.12E+01  4.10E+01  7.13E-03  4.06E+01
 
 TH 5
+       -1.86E+01 -3.61E+01 -7.98E-03 -8.55E+00  1.14E+01
 
 TH 6
+        2.80E+01 -1.13E+01  3.87E-03  9.01E+00  2.20E+00  4.59E+01
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -1.74E+00 -1.97E+01  6.75E-03  4.90E+00  2.89E+00  2.80E+00  0.00E+00  0.00E+00  7.08E+00
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -2.79E+01 -2.03E+01  5.58E-03 -3.78E+00  5.64E+00 -5.39E+00  0.00E+00  0.00E+00  7.72E+00  0.00E+00  1.86E+02
 
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
 
 Elapsed finaloutput time in seconds:     0.03
 #CPUT: Total CPU Time in Seconds,       74.387
Stop Time:
Thu Sep 30 09:13:29 CDT 2021
