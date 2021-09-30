Wed Sep 29 11:56:32 CDT 2021
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
$DATA ../../../../data/spa/A1/dat15.csv ignore=@
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1423.00592051410        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.3197E+02 -4.9139E+01  2.5000E+00 -6.1324E+01  1.0130E+02  5.5945E+01 -2.5327E+01 -1.2634E+01 -5.6760E+01 -3.4601E+01
            -4.0833E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1529.25466963527        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.0823E+00  1.0179E+00  9.4777E-01  1.0726E+00  9.3192E-01  8.8968E-01  1.0653E+00  1.0100E+00  1.1643E+00  1.0038E+00
             1.9180E+00
 PARAMETER:  1.7908E-01  1.1774E-01  4.6360E-02  1.7009E-01  2.9492E-02 -1.6889E-02  1.6321E-01  1.0995E-01  2.5212E-01  1.0381E-01
             7.5128E-01
 GRADIENT:   2.8126E+02 -9.9507E+00 -1.6566E+01  2.0141E+01  2.1298E+01 -1.0652E+01  1.5072E+00  4.2213E+00  1.1828E+01  5.9344E+00
             4.9891E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1530.09352200265        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      266
 NPARAMETR:  1.0756E+00  8.1847E-01  1.0703E+00  1.2612E+00  8.8866E-01  9.1900E-01  1.4240E+00  8.6515E-01  1.0330E+00  8.9951E-01
             1.9050E+00
 PARAMETER:  1.7289E-01 -1.0032E-01  1.6791E-01  3.3204E-01 -1.8038E-02  1.5528E-02  4.5349E-01 -4.4856E-02  1.3250E-01 -5.9060E-03
             7.4446E-01
 GRADIENT:   5.4878E+01  2.2554E+01  2.5064E+00  4.8928E+01  6.9238E+00 -5.7507E+00  4.3390E+00 -2.1052E+00  2.3505E+00 -6.9179E+00
             3.2475E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1534.42000048507        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      443
 NPARAMETR:  1.0535E+00  8.2268E-01  8.8383E-01  1.1954E+00  7.9195E-01  9.3407E-01  1.3681E+00  7.3510E-01  1.0338E+00  9.1711E-01
             1.7097E+00
 PARAMETER:  1.5211E-01 -9.5192E-02 -2.3494E-02  2.7850E-01 -1.3326E-01  3.1795E-02  4.1340E-01 -2.0775E-01  1.3321E-01  1.3470E-02
             6.3634E-01
 GRADIENT:   4.3555E+00  5.8007E+00  7.2467E+00 -3.0134E+00 -1.2102E+01 -8.3560E-02  5.0265E-01  3.8014E-02  6.7625E-02  6.7841E-02
             2.9427E-02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1535.23315166014        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      619
 NPARAMETR:  1.0459E+00  5.0786E-01  1.0407E+00  1.3958E+00  7.6811E-01  9.3095E-01  1.6124E+00  8.6751E-01  9.6772E-01  9.7215E-01
             1.7196E+00
 PARAMETER:  1.4491E-01 -5.7755E-01  1.3991E-01  4.3349E-01 -1.6382E-01  2.8454E-02  5.7772E-01 -4.2123E-02  6.7189E-02  7.1755E-02
             6.4208E-01
 GRADIENT:  -1.7840E+00  3.9852E+00  7.1585E+00 -3.3929E-01 -1.2481E+01  3.5776E-01  1.3217E-01  2.7822E-01  2.6260E+00  2.1369E+00
             3.7730E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1536.08398471228        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      796
 NPARAMETR:  1.0419E+00  2.2596E-01  1.0340E+00  1.5670E+00  6.9886E-01  9.2632E-01  2.0873E+00  8.4931E-01  8.8917E-01  9.6801E-01
             1.6963E+00
 PARAMETER:  1.4108E-01 -1.3874E+00  1.3339E-01  5.4914E-01 -2.5831E-01  2.3469E-02  8.3587E-01 -6.3328E-02 -1.7462E-02  6.7482E-02
             6.2846E-01
 GRADIENT:   1.4889E+00  3.2100E+00  4.1718E+00  2.0270E+01 -8.6874E+00 -2.7222E-01 -7.8456E-01 -2.5880E-01 -1.5977E+00  8.5025E-01
            -2.1276E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1536.56913697845        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      975
 NPARAMETR:  1.0365E+00  8.3229E-02  1.0589E+00  1.6461E+00  6.8192E-01  9.2403E-01  3.1690E+00  9.0632E-01  8.5472E-01  9.5375E-01
             1.6974E+00
 PARAMETER:  1.3585E-01 -2.3862E+00  1.5725E-01  5.9843E-01 -2.8284E-01  2.0990E-02  1.2534E+00  1.6357E-03 -5.6983E-02  5.2644E-02
             6.2912E-01
 GRADIENT:  -2.9920E+00  9.2877E-01  1.0367E+00  1.5879E+01 -2.4449E+00 -2.4727E-01 -3.1555E-01 -6.4418E-02 -1.9744E+00 -7.1575E-01
            -1.5681E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1536.93008809122        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1155
 NPARAMETR:  1.0359E+00  1.4336E-02  9.8784E-01  1.6658E+00  6.3709E-01  9.2361E-01  7.5502E+00  8.0636E-01  8.4531E-01  9.5843E-01
             1.7002E+00
 PARAMETER:  1.3530E-01 -4.1450E+00  8.7768E-02  6.1030E-01 -3.5084E-01  2.0535E-02  2.1216E+00 -1.1522E-01 -6.8047E-02  5.7537E-02
             6.3074E-01
 GRADIENT:  -8.3106E-02  1.1825E-01  1.8178E+00  4.3299E+00 -4.3619E+00  7.2453E-02 -5.3243E-02  4.7368E-02 -1.3477E-01  6.0959E-01
             5.4311E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1536.97430343966        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:     1316
 NPARAMETR:  1.0353E+00  1.0000E-02  9.9201E-01  1.6643E+00  6.3956E-01  9.2351E-01  9.9635E+00  8.0910E-01  8.4568E-01  9.5747E-01
             1.7001E+00
 PARAMETER:  1.3471E-01 -4.5466E+00  9.1974E-02  6.0938E-01 -3.4697E-01  2.0424E-02  2.3989E+00 -1.1183E-01 -6.7617E-02  5.6544E-02
             6.3070E-01
 GRADIENT:   1.8662E+02  6.3480E-03  1.9652E+00  3.7857E+02  1.3026E+01  1.0460E+01  1.5904E-01  2.4197E-02  6.2610E+00  6.3431E-01
             5.5736E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1536.98244500940        NO. OF FUNC. EVALS.: 164
 CUMULATIVE NO. OF FUNC. EVALS.:     1480
 NPARAMETR:  1.0357E+00  1.0000E-02  9.9143E-01  1.6624E+00  6.3962E-01  9.2328E-01  1.2141E+01  8.1182E-01  8.4417E-01  9.5621E-01
             1.6985E+00
 PARAMETER:  1.3509E-01 -4.5534E+00  9.1397E-02  6.0824E-01 -3.4688E-01  2.0173E-02  2.5966E+00 -1.0848E-01 -6.9396E-02  5.5220E-02
             6.2977E-01
 GRADIENT:   6.9024E-02  0.0000E+00  1.3868E-01 -6.9442E+00  1.9498E-01  1.4345E-02  2.0550E-02 -1.3673E-02  3.5540E-02  5.9309E-02
             1.7179E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1536.98335938533        NO. OF FUNC. EVALS.: 193
 CUMULATIVE NO. OF FUNC. EVALS.:     1673             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0360E+00  1.0000E-02  9.9024E-01  1.6617E+00  6.3900E-01  9.2330E-01  1.1918E+01  8.1201E-01  8.4429E-01  9.5488E-01
             1.6976E+00
 PARAMETER:  1.3537E-01 -4.5534E+00  9.0189E-02  6.0787E-01 -3.4784E-01  2.0199E-02  2.5781E+00 -1.0824E-01 -6.9257E-02  5.3826E-02
             6.2924E-01
 GRADIENT:   1.8992E+02  0.0000E+00  1.6427E+00  3.7493E+02  1.4022E+01  1.0387E+01  3.0150E-01  8.1527E-02  5.7423E+00  4.3923E-01
             4.9761E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1536.98359032358        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     1861
 NPARAMETR:  1.0360E+00  1.0000E-02  9.8957E-01  1.6616E+00  6.3872E-01  9.2330E-01  1.1794E+01  8.1069E-01  8.4422E-01  9.5493E-01
             1.6977E+00
 PARAMETER:  1.3540E-01 -4.5534E+00  8.9511E-02  6.0779E-01 -3.4830E-01  2.0197E-02  2.5676E+00 -1.0987E-01 -6.9340E-02  5.3881E-02
             6.2928E-01
 GRADIENT:   8.9804E-01  0.0000E+00  2.1872E-01 -7.6661E+00  1.5932E-01  1.9602E-02  4.8217E-03 -2.0520E-02  6.7792E-04 -4.9426E-02
            -5.3868E-02

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1536.98378564937        NO. OF FUNC. EVALS.: 192
 CUMULATIVE NO. OF FUNC. EVALS.:     2053             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0360E+00  1.0000E-02  9.8872E-01  1.6615E+00  6.3841E-01  9.2331E-01  1.1840E+01  8.1149E-01  8.4428E-01  9.5540E-01
             1.6979E+00
 PARAMETER:  1.3539E-01 -4.5534E+00  8.8656E-02  6.0772E-01 -3.4877E-01  2.0204E-02  2.5715E+00 -1.0888E-01 -6.9273E-02  5.4380E-02
             6.2938E-01
 GRADIENT:   1.8989E+02  0.0000E+00  1.4318E+00  3.7460E+02  1.4158E+01  1.0385E+01  2.9291E-01  1.5280E-01  5.7155E+00  5.7653E-01
             5.1279E+00

0ITERATION NO.:   64    OBJECTIVE VALUE:  -1536.98388507622        NO. OF FUNC. EVALS.: 134
 CUMULATIVE NO. OF FUNC. EVALS.:     2187
 NPARAMETR:  1.0360E+00  1.0000E-02  9.8872E-01  1.6615E+00  6.3834E-01  9.2330E-01  1.1838E+01  8.0985E-01  8.4428E-01  9.5495E-01
             1.6977E+00
 PARAMETER:  1.3539E-01 -4.5534E+00  8.8652E-02  6.0771E-01 -3.4888E-01  2.0201E-02  2.5713E+00 -1.1090E-01 -6.9268E-02  5.3903E-02
             6.2925E-01
 GRADIENT:   8.8125E-01  0.0000E+00  1.9243E-01 -7.6608E+00  1.4229E-01  2.0011E-02  6.4407E-03 -9.0369E-03  1.3855E-02 -2.6107E-02
            -4.0341E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     2187
 NO. OF SIG. DIGITS IN FINAL EST.:  2.2
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.3558E-04  1.3311E-04 -1.1774E-02 -7.4320E-03 -2.4031E-02
 SE:             2.9524E-02  1.6853E-03  1.3488E-02  2.8489E-02  2.1427E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9634E-01  9.3704E-01  3.8272E-01  7.9419E-01  2.6206E-01

 ETASHRINKSD(%)  1.0905E+00  9.4354E+01  5.4812E+01  4.5575E+00  2.8218E+01
 ETASHRINKVR(%)  2.1692E+00  9.9681E+01  7.9581E+01  8.9073E+00  4.8473E+01
 EBVSHRINKSD(%)  1.3176E+00  9.4593E+01  5.5726E+01  4.6036E+00  2.6959E+01
 EBVSHRINKVR(%)  2.6179E+00  9.9708E+01  8.0398E+01  8.9953E+00  4.6650E+01
 RELATIVEINF(%)  8.7566E+01  1.0891E-02  1.9830E+00  5.2281E+00  2.7923E+00
 EPSSHRINKSD(%)  3.9870E+01
 EPSSHRINKVR(%)  6.3844E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1536.9838850762201     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -801.83305851248190     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    28.23
0S MATRIX ALGORITHMICALLY SINGULAR
0S MATRIX IS OUTPUT
0INVERSE COVARIANCE MATRIX SET TO RS*RMAT, WHERE S* IS A PSEUDO INVERSE OF S
 Elapsed covariance  time in seconds:     5.80
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1536.984       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.04E+00  1.00E-02  9.89E-01  1.66E+00  6.38E-01  9.23E-01  1.18E+01  8.10E-01  8.44E-01  9.55E-01  1.70E+00
 


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
 
         2.97E-02  0.00E+00  5.52E-01  1.11E-01  2.47E-01  7.79E-02  2.00E+00  8.69E-01  1.05E-01  1.19E-01  1.50E-01
 


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
+        8.82E-04
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -6.45E-04  0.00E+00  3.05E-01
 
 TH 4
+       -1.01E-04  0.00E+00  5.68E-02  1.24E-02
 
 TH 5
+       -2.78E-04  0.00E+00  1.35E-01  2.52E-02  6.08E-02
 
 TH 6
+       -7.42E-04  0.00E+00 -3.86E-04 -1.52E-04 -1.18E-04  6.07E-03
 
 TH 7
+        5.67E-03  0.00E+00 -7.87E-01 -1.72E-01 -3.51E-01 -9.78E-03  3.99E+00
 
 TH 8
+       -2.32E-03  0.00E+00  4.52E-01  8.54E-02  2.01E-01 -1.35E-03 -1.25E+00  7.56E-01
 
 TH 9
+       -8.76E-05  0.00E+00 -2.74E-02 -4.42E-03 -1.21E-02 -4.76E-04  1.04E-01 -4.16E-02  1.10E-02
 
 TH10
+        3.41E-04  0.00E+00 -9.91E-03 -2.39E-03 -3.68E-03  4.17E-04  8.64E-02 -3.28E-02  8.59E-04  1.42E-02
 
 TH11
+        1.69E-03  0.00E+00 -2.96E-02 -6.18E-03 -1.25E-02  1.98E-04  1.02E-01 -6.01E-02  1.78E-03  5.22E-03  2.26E-02
 
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
+        2.97E-02
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -3.93E-02  0.00E+00  5.52E-01
 
 TH 4
+       -3.05E-02  0.00E+00  9.25E-01  1.11E-01
 
 TH 5
+       -3.80E-02  0.00E+00  9.95E-01  9.19E-01  2.47E-01
 
 TH 6
+       -3.21E-01  0.00E+00 -8.96E-03 -1.76E-02 -6.12E-03  7.79E-02
 
 TH 7
+        9.55E-02  0.00E+00 -7.13E-01 -7.72E-01 -7.12E-01 -6.28E-02  2.00E+00
 
 TH 8
+       -9.01E-02  0.00E+00  9.42E-01  8.83E-01  9.36E-01 -2.00E-02 -7.22E-01  8.69E-01
 
 TH 9
+       -2.81E-02  0.00E+00 -4.72E-01 -3.79E-01 -4.67E-01 -5.83E-02  4.95E-01 -4.56E-01  1.05E-01
 
 TH10
+        9.63E-02  0.00E+00 -1.50E-01 -1.80E-01 -1.25E-01  4.49E-02  3.62E-01 -3.16E-01  6.85E-02  1.19E-01
 
 TH11
+        3.79E-01  0.00E+00 -3.57E-01 -3.70E-01 -3.37E-01  1.69E-02  3.40E-01 -4.60E-01  1.13E-01  2.91E-01  1.50E-01
 
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
+        1.23E+03
 
 TH 2
+       -1.10E-13 -1.87E-28
 
 TH 3
+       -1.40E+01 -2.17E-14  3.21E+02
 
 TH 4
+        4.17E+01  1.68E-13 -6.16E+01  4.96E+02
 
 TH 5
+        1.09E+00 -3.65E-14 -7.04E+02 -3.99E+01  1.61E+03
 
 TH 6
+        1.15E+02  2.37E-13 -5.33E+00  3.24E+01 -1.78E+01  1.68E+02
 
 TH 7
+        9.10E-04 -7.01E-19  2.67E-04 -5.15E-03  1.08E-03  8.54E-04  1.09E-06
 
 TH 8
+        5.75E+00  6.13E-15 -3.58E+00 -2.57E+00  8.45E+00  5.02E+00  2.41E-04  2.57E-01
 
 TH 9
+       -7.44E+00 -7.83E-15 -8.48E+00  2.76E+00  1.80E+01 -2.03E+00  9.77E-03  1.64E+00  9.40E+01
 
 TH10
+        8.83E+00  1.21E-14  2.72E+01 -4.92E+00 -6.08E+01  9.89E+00  6.44E-04  1.14E-01  4.37E+00  3.24E+00
 
 TH11
+        6.99E+00  1.85E-14  2.97E-02 -1.56E+01  3.54E+00  1.71E+01  1.61E-03  8.80E-01  1.23E+01  1.81E+00  4.18E+00
 
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
+        1.18E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+        4.12E+00  0.00E+00  3.43E+02
 
 TH 4
+       -2.47E+01  0.00E+00  6.42E+00  5.14E+02
 
 TH 5
+        5.35E+00  0.00E+00 -6.51E+02 -1.07E+02  1.62E+03
 
 TH 6
+       -1.89E+02  0.00E+00 -3.01E+00 -1.89E+01  1.44E+01  2.98E+02
 
 TH 7
+        1.21E-02  0.00E+00  2.11E-03 -2.96E-02 -1.40E-02 -5.66E-03  1.58E-05
 
 TH 8
+       -1.98E+00  0.00E+00 -3.93E+01 -5.21E+00  2.38E+01 -4.52E-01  1.20E-03  1.98E+01
 
 TH 9
+       -1.59E+00  0.00E+00  7.16E+00  1.23E+01 -7.26E+01 -3.09E+01  6.83E-02 -1.09E+01  5.18E+02
 
 TH10
+        1.55E+01  0.00E+00 -1.91E+01 -3.55E+01 -6.71E+01  3.78E+00  1.47E-02  1.95E+01  2.98E+01  6.24E+01
 
 TH11
+        1.03E+02  0.00E+00 -4.52E+01 -7.89E+01  5.78E+00  1.25E+01  1.97E-02  1.06E+01  7.28E+01  5.13E+01  1.46E+02
 
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
 
 Elapsed finaloutput time in seconds:     0.04
1THERE ARE ERROR MESSAGES IN FILE PRDERR                                                                  
 #CPUT: Total CPU Time in Seconds,       34.093
Stop Time:
Wed Sep 29 11:57:08 CDT 2021
