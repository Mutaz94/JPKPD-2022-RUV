Wed Sep 29 12:25:01 CDT 2021
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
$DATA ../../../../data/spa/A1/dat86.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m86.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1543.83062802790        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.0382E+02  6.1081E+00 -3.2827E+01  8.1519E+01  1.3656E+02  5.7572E+01  3.8791E+00 -6.0202E+00 -2.4755E+00 -1.5605E+01
            -2.5534E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1595.38544052474        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  1.0252E+00  9.8828E-01  9.5531E-01  1.0007E+00  9.1197E-01  9.6826E-01  9.6886E-01  1.0216E+00  1.0124E+00  9.3112E-01
             1.4715E+00
 PARAMETER:  1.2487E-01  8.8211E-02  5.4283E-02  1.0073E-01  7.8566E-03  6.7742E-02  6.8363E-02  1.2140E-01  1.1233E-01  2.8633E-02
             4.8630E-01
 GRADIENT:   2.5265E+02 -7.2349E+00 -1.2029E+01  1.7339E+01  4.3275E+01  2.0114E+01  4.5517E+00 -3.1123E-01  3.4234E+00  8.8737E+00
             1.9132E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1597.78046503652        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  1.0164E+00  1.0561E+00  6.9507E-01  9.6078E-01  7.8240E-01  9.5834E-01  1.0372E+00  8.5357E-01  9.8018E-01  7.1366E-01
             1.4508E+00
 PARAMETER:  1.1623E-01  1.5462E-01 -2.6374E-01  5.9990E-02 -1.4539E-01  5.7443E-02  1.3652E-01 -5.8325E-02  7.9983E-02 -2.3735E-01
             4.7213E-01
 GRADIENT:   2.0906E+02  4.3638E+01  2.7346E-01  4.7491E+01  1.4654E+01  1.4782E+01  9.4323E+00  1.4766E+00 -2.6581E+00  5.0834E+00
            -4.8140E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1597.98299336905        NO. OF FUNC. EVALS.: 125
 CUMULATIVE NO. OF FUNC. EVALS.:      281
 NPARAMETR:  1.0160E+00  1.0786E+00  6.9285E-01  9.4829E-01  7.8742E-01  9.6116E-01  1.0164E+00  8.3396E-01  1.0345E+00  7.1030E-01
             1.4518E+00
 PARAMETER:  1.1586E-01  1.7568E-01 -2.6694E-01  4.6907E-02 -1.3899E-01  6.0384E-02  1.1631E-01 -8.1564E-02  1.3394E-01 -2.4207E-01
             4.7283E-01
 GRADIENT:  -4.2761E+01  1.4696E+01  1.5050E+00  1.2962E+01 -1.3250E+00 -6.8519E+00  7.5884E+00  7.7071E-01  1.3000E+00  3.6141E+00
            -8.1292E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1599.89990746511        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      457
 NPARAMETR:  1.0361E+00  9.7854E-01  4.6323E-01  9.5413E-01  6.1477E-01  9.7676E-01  1.0681E+00  4.2462E-01  9.5363E-01  5.0449E-01
             1.4707E+00
 PARAMETER:  1.3548E-01  7.8309E-02 -6.6954E-01  5.3046E-02 -3.8651E-01  7.6481E-02  1.6589E-01 -7.5657E-01  5.2515E-02 -5.8420E-01
             4.8576E-01
 GRADIENT:   2.6565E-01  1.2131E+00  9.3127E-01  4.0249E+00 -2.5954E+00 -6.1112E-01  1.7307E+00 -1.1479E-01  1.2095E+00  1.0919E+00
             5.8103E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1601.18977618667        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      635
 NPARAMETR:  1.0320E+00  1.1872E+00  2.9121E-01  7.8026E-01  5.9175E-01  9.7804E-01  9.0047E-01  2.9308E-01  1.0027E+00  3.2545E-01
             1.4713E+00
 PARAMETER:  1.3152E-01  2.7162E-01 -1.1337E+00 -1.4813E-01 -4.2467E-01  7.7792E-02 -4.8346E-03 -1.1273E+00  1.0270E-01 -1.0226E+00
             4.8615E-01
 GRADIENT:   3.2799E+00  6.9034E+00  5.9758E+00 -2.7542E+00 -9.2627E+00  4.8818E-01 -8.1071E-01 -1.4036E-01 -5.0354E+00 -2.5484E-01
            -6.2870E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1601.74944254545        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      810
 NPARAMETR:  1.0236E+00  1.3505E+00  2.2706E-01  6.6588E-01  6.2748E-01  9.7275E-01  8.2183E-01  2.4610E-01  1.1242E+00  2.5685E-01
             1.4785E+00
 PARAMETER:  1.2331E-01  4.0045E-01 -1.3825E+00 -3.0664E-01 -3.6604E-01  7.2369E-02 -9.6219E-02 -1.3020E+00  2.1709E-01 -1.2593E+00
             4.9100E-01
 GRADIENT:  -4.6551E+00 -2.1560E+00 -1.5442E-01  1.4670E+00  2.8620E+00  2.9764E-02  4.5606E-01 -1.0598E-01 -6.2082E-02 -5.7124E-01
            -1.5746E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1601.95636241809        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      989
 NPARAMETR:  1.0263E+00  1.6041E+00  1.8412E-01  5.2560E-01  7.2839E-01  9.6820E-01  7.2452E-01  1.8884E-01  1.3218E+00  3.3818E-01
             1.5088E+00
 PARAMETER:  1.2593E-01  5.7254E-01 -1.5921E+00 -5.4322E-01 -2.1692E-01  6.7680E-02 -2.2224E-01 -1.5669E+00  3.7897E-01 -9.8419E-01
             5.1133E-01
 GRADIENT:  -2.1339E+01  1.6757E+01  2.8080E+00  4.4577E+00 -3.7850E+00 -3.5480E+00 -2.7041E+00 -2.3666E-01 -1.1765E+00 -1.8058E+00
             1.4138E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1602.27801445841        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1166
 NPARAMETR:  1.0286E+00  1.6984E+00  1.6675E-01  4.6676E-01  7.7591E-01  9.7156E-01  6.9865E-01  1.5594E-01  1.4268E+00  4.6722E-01
             1.4873E+00
 PARAMETER:  1.2815E-01  6.2966E-01 -1.6913E+00 -6.6194E-01 -1.5372E-01  7.1146E-02 -2.5861E-01 -1.7583E+00  4.5542E-01 -6.6094E-01
             4.9695E-01
 GRADIENT:   9.3665E+00 -4.1993E+00  6.2933E-01 -1.9662E+00 -1.7981E+00  1.4116E+00  5.7721E-01 -5.3044E-02 -9.5131E-01  5.8597E-01
             9.7344E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1602.28778858340        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1342
 NPARAMETR:  1.0283E+00  1.7605E+00  1.5395E-01  4.3419E-01  8.0171E-01  9.7127E-01  6.8428E-01  1.4371E-01  1.5043E+00  4.9182E-01
             1.4813E+00
 PARAMETER:  1.2788E-01  6.6560E-01 -1.7711E+00 -7.3426E-01 -1.2100E-01  7.0850E-02 -2.7938E-01 -1.8399E+00  5.0836E-01 -6.0964E-01
             4.9294E-01
 GRADIENT:   9.0210E+00  6.6815E+00  1.9002E-01  1.7860E+00 -5.4749E+00  1.4555E+00  7.6699E-02 -5.9378E-02 -2.3861E-01  3.5522E-01
            -1.2050E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1602.28843513422        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1521
 NPARAMETR:  1.0280E+00  1.7785E+00  1.5056E-01  4.2405E-01  8.1084E-01  9.7073E-01  6.8056E-01  1.4145E-01  1.5293E+00  5.0007E-01
             1.4815E+00
 PARAMETER:  1.2758E-01  6.7578E-01 -1.7934E+00 -7.5789E-01 -1.0969E-01  7.0289E-02 -2.8484E-01 -1.8558E+00  5.2478E-01 -5.9300E-01
             4.9306E-01
 GRADIENT:   8.1918E+00  6.6030E+00 -3.6338E-02  2.2572E+00 -4.5463E+00  1.2898E+00  3.6691E-02 -6.1552E-02 -1.4443E-01  2.9390E-01
            -1.2654E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1602.30917729734        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1705
 NPARAMETR:  1.0227E+00  1.7781E+00  1.4959E-01  4.1920E-01  8.1631E-01  9.6503E-01  6.7952E-01  1.4626E-01  1.5381E+00  4.9592E-01
             1.4861E+00
 PARAMETER:  1.2242E-01  6.7556E-01 -1.7998E+00 -7.6942E-01 -1.0297E-01  6.4400E-02 -2.8636E-01 -1.8224E+00  5.3055E-01 -6.0135E-01
             4.9618E-01
 GRADIENT:  -3.8474E+00 -1.1722E+01 -1.0714E+00  4.5897E-01  7.0622E+00 -8.9468E-01 -1.0257E-01 -6.6544E-02 -3.9855E-01 -1.1589E-01
            -1.8736E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1602.50984399270        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1886
 NPARAMETR:  1.0256E+00  1.7924E+00  1.4755E-01  4.1275E-01  8.2010E-01  9.6836E-01  6.7727E-01  3.5457E-01  1.5571E+00  5.0401E-01
             1.4785E+00
 PARAMETER:  1.2530E-01  6.8357E-01 -1.8136E+00 -7.8491E-01 -9.8327E-02  6.7850E-02 -2.8969E-01 -9.3685E-01  5.4285E-01 -5.8516E-01
             4.9104E-01
 GRADIENT:   3.1392E+00 -2.5173E+00  6.1275E-01  1.7403E-01 -1.1790E+00  5.0341E-01  7.4668E-02 -5.1670E-01  3.9660E-01  2.1120E-01
            -8.4874E-01

0ITERATION NO.:   63    OBJECTIVE VALUE:  -1602.61095689077        NO. OF FUNC. EVALS.: 100
 CUMULATIVE NO. OF FUNC. EVALS.:     1986
 NPARAMETR:  1.0259E+00  1.7899E+00  1.4779E-01  4.1113E-01  8.2031E-01  9.6760E-01  6.7661E-01  4.1583E-01  1.5579E+00  5.0193E-01
             1.4759E+00
 PARAMETER:  1.2534E-01  6.8343E-01 -1.8153E+00 -7.8739E-01 -9.8581E-02  6.7738E-02 -2.9089E-01 -7.7605E-01  5.4434E-01 -5.8345E-01
             4.8838E-01
 GRADIENT:  -9.5047E+04  1.7432E+04 -6.5303E+03  1.5114E+04 -1.5119E+00  3.4932E-01 -1.4623E-01  1.5336E+04  1.0922E+04  2.7504E-01
            -1.2288E+04
 NUMSIGDIG:         2.3         2.3         2.3         2.3         1.9         1.7         2.7         2.3         2.3         1.6
                    2.3

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1986
 NO. OF SIG. DIGITS IN FINAL EST.:  1.6

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -6.0929E-04 -2.1345E-02 -3.5530E-03  2.1159E-02 -2.6501E-02
 SE:             2.9661E-02  2.5843E-02  3.9264E-03  2.4002E-02  1.5022E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8361E-01  4.0883E-01  3.6552E-01  3.7802E-01  7.7703E-02

 ETASHRINKSD(%)  6.3159E-01  1.3422E+01  8.6846E+01  1.9591E+01  4.9675E+01
 ETASHRINKVR(%)  1.2592E+00  2.5043E+01  9.8270E+01  3.5344E+01  7.4674E+01
 EBVSHRINKSD(%)  9.0143E-01  1.4212E+01  8.8887E+01  1.8240E+01  5.0039E+01
 EBVSHRINKVR(%)  1.7947E+00  2.6404E+01  9.8765E+01  3.3154E+01  7.5039E+01
 RELATIVEINF(%)  9.6847E+01  9.4823E+00  3.6461E-01  8.6723E+00  3.5575E+00
 EPSSHRINKSD(%)  4.1064E+01
 EPSSHRINKVR(%)  6.5265E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1602.6109568907668     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -867.46013032702865     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    23.65
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.94
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1602.611       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  1.79E+00  1.47E-01  4.12E-01  8.20E-01  9.68E-01  6.76E-01  4.16E-01  1.56E+00  5.05E-01  1.47E+00
 


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
+        1.80E+07
 
 TH 2
+        1.99E+02  1.99E+05
 
 TH 3
+       -1.10E+03  6.43E+02  4.13E+06
 
 TH 4
+        8.18E+02 -1.44E+03  5.96E+03  2.83E+06
 
 TH 5
+        2.83E+07 -2.07E+02 -3.31E+03  2.22E+03  1.63E+03
 
 TH 6
+       -1.55E+03  1.59E+02 -7.52E+02  6.06E+02 -5.45E+00  2.03E+02
 
 TH 7
+       -6.40E+02  7.32E+01 -3.58E+02  2.61E+02  7.54E+00  1.51E+00  2.49E+02
 
 TH 8
+        8.32E+02 -1.48E+03  6.69E+03 -6.37E+03  1.43E+03  6.14E+02  2.57E+02  2.85E+06
 
 TH 9
+        3.24E+02 -1.16E+03  5.15E+03  1.08E+06  5.15E+02  2.35E+02  1.10E+02  1.08E+06  4.12E+05
 
 TH10
+        7.86E+06  2.12E+00 -1.08E+02  4.67E+01 -8.20E+01  1.36E+00  2.22E+01  6.97E+01  3.44E+01  6.79E+01
 
 TH11
+       -3.87E+02 -4.98E+03  2.27E+04  2.86E+03 -6.75E+02 -2.74E+02 -1.01E+02  2.53E+03  1.97E+03 -4.14E+00  5.83E+05
 
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
 #CPUT: Total CPU Time in Seconds,       29.634
Stop Time:
Wed Sep 29 12:25:33 CDT 2021
