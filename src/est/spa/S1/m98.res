Wed Sep 29 14:45:00 CDT 2021
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
$DATA ../../../../data/spa/S1/dat98.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m98.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1681.59732229113        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.6539E+02 -1.7029E+01 -3.5962E+01  5.3628E+01  6.1370E+01  6.1948E+01 -4.3231E+00  2.6103E+00  2.0077E+01  4.2337E+00
             4.9130E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1684.90250277546        NO. OF FUNC. EVALS.: 119
 CUMULATIVE NO. OF FUNC. EVALS.:      132
 NPARAMETR:  9.6832E-01  1.0790E+00  1.0442E+00  1.0071E+00  9.4556E-01  8.8350E-01  1.0175E+00  9.9592E-01  9.6226E-01  9.8174E-01
             8.9520E-01
 PARAMETER:  6.7812E-02  1.7607E-01  1.4323E-01  1.0703E-01  4.4024E-02 -2.3869E-02  1.1737E-01  9.5912E-02  6.1533E-02  8.1574E-02
            -1.0712E-02
 GRADIENT:  -1.2241E+01  6.1592E+01  3.3966E+01  4.7346E+01 -8.0909E+01 -2.5954E+01 -8.2379E+00 -5.1643E+00  3.8546E+00  2.8931E+00
             5.3870E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1687.72676936880        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      308
 NPARAMETR:  9.7230E-01  1.0851E+00  1.2307E+00  1.0017E+00  1.0356E+00  9.2748E-01  1.2632E+00  1.4592E+00  9.1209E-01  9.2075E-01
             8.6086E-01
 PARAMETER:  7.1908E-02  1.8163E-01  3.0760E-01  1.0171E-01  1.3501E-01  2.4718E-02  3.3367E-01  4.7786E-01  7.9863E-03  1.7436E-02
            -4.9825E-02
 GRADIENT:   2.4439E+00  4.1577E+01  1.5225E+01  3.0958E+01 -2.6524E+01 -5.0070E+00  1.2434E+01 -7.7185E-01  5.5089E+00 -5.4845E+00
            -8.4400E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1690.35508287808        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      488
 NPARAMETR:  9.7186E-01  1.2213E+00  8.8492E-01  8.7604E-01  1.0034E+00  9.4100E-01  1.0581E+00  1.0849E+00  9.5409E-01  9.3011E-01
             8.7576E-01
 PARAMETER:  7.1454E-02  2.9995E-01 -2.2255E-02 -3.2339E-02  1.0340E-01  3.9187E-02  1.5643E-01  1.8149E-01  5.3008E-02  2.7543E-02
            -3.2669E-02
 GRADIENT:  -4.2652E+00  1.7794E+00 -3.4115E+00  8.1032E+00  3.4582E+00 -5.5333E-02  1.1408E+00  6.4780E-01  7.1982E-01  2.5578E-01
            -4.7413E-02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1690.58052942319        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      663
 NPARAMETR:  9.7565E-01  1.4696E+00  7.1640E-01  7.1154E-01  1.0500E+00  9.4216E-01  9.2230E-01  1.0203E+00  1.0815E+00  9.4940E-01
             8.7752E-01
 PARAMETER:  7.5345E-02  4.8499E-01 -2.3351E-01 -2.4033E-01  1.4879E-01  4.0420E-02  1.9116E-02  1.2014E-01  1.7831E-01  4.8074E-02
            -3.0660E-02
 GRADIENT:   2.7428E+00  5.6483E+00  3.3889E+00  9.8084E-01 -5.2537E+00  8.3253E-02 -5.7924E-01 -7.1654E-01 -7.5075E-01 -2.9733E-01
            -9.8153E-02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1690.62755556188        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      845
 NPARAMETR:  9.7548E-01  1.5306E+00  6.5776E-01  6.6753E-01  1.0616E+00  9.4242E-01  8.9831E-01  1.0209E+00  1.1239E+00  9.5186E-01
             8.7736E-01
 PARAMETER:  7.5179E-02  5.2569E-01 -3.1892E-01 -3.0416E-01  1.5981E-01  4.0691E-02 -7.2407E-03  1.2068E-01  2.1679E-01  5.0667E-02
            -3.0837E-02
 GRADIENT:   1.6873E+00 -2.7169E+00 -1.1181E-01  4.5102E-01  6.5605E-01  9.3306E-02  2.7647E-02 -1.6116E-02  1.2395E-01  3.2146E-01
             2.2512E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1690.62882472833        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1026
 NPARAMETR:  9.7548E-01  1.5308E+00  6.5740E-01  6.6680E-01  1.0611E+00  9.4254E-01  8.9815E-01  1.0260E+00  1.1233E+00  9.4854E-01
             8.7690E-01
 PARAMETER:  7.5173E-02  5.2582E-01 -3.1946E-01 -3.0527E-01  1.5929E-01  4.0819E-02 -7.4148E-03  1.2572E-01  2.1626E-01  4.7167E-02
            -3.1362E-02
 GRADIENT:   1.6894E+00 -3.1891E+00  2.6292E-01 -4.7641E-01  2.2095E-01  1.4528E-01 -1.5487E-02 -3.2756E-02 -2.9136E-02 -2.8083E-02
            -1.8247E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1690.62902374028        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1207
 NPARAMETR:  9.7560E-01  1.5308E+00  6.5676E-01  6.6672E-01  1.0605E+00  9.4260E-01  8.9805E-01  1.0233E+00  1.1230E+00  9.4777E-01
             8.7683E-01
 PARAMETER:  7.5302E-02  5.2578E-01 -3.2044E-01 -3.0539E-01  1.5876E-01  4.0891E-02 -7.5268E-03  1.2301E-01  2.1601E-01  4.6357E-02
            -3.1443E-02
 GRADIENT:   2.0042E+00 -3.1463E+00  4.0959E-01 -6.6842E-01 -1.5686E-02  1.7191E-01 -6.6827E-02 -6.0897E-02 -7.5791E-02 -8.8852E-02
            -6.8198E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1690.62941075791        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     1392
 NPARAMETR:  9.7561E-01  1.5308E+00  6.5607E-01  6.6671E-01  1.0602E+00  9.4261E-01  8.9815E-01  1.0233E+00  1.1230E+00  9.4770E-01
             8.7685E-01
 PARAMETER:  7.5304E-02  5.2582E-01 -3.2149E-01 -3.0539E-01  1.5846E-01  4.0895E-02 -7.4218E-03  1.2307E-01  2.1604E-01  4.6280E-02
            -3.1424E-02
 GRADIENT:   1.9963E+00 -3.1396E+00  2.9492E-01 -5.0645E-01  3.0266E-02  1.7064E-01 -3.4350E-02 -3.1107E-02 -2.8106E-02 -2.2368E-02
            -3.6843E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1690.62967623471        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     1580             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7561E-01  1.5306E+00  6.5488E-01  6.6685E-01  1.0598E+00  9.4261E-01  8.9840E-01  1.0252E+00  1.1232E+00  9.4749E-01
             8.7690E-01
 PARAMETER:  7.5306E-02  5.2566E-01 -3.2330E-01 -3.0519E-01  1.5810E-01  4.0900E-02 -7.1354E-03  1.2490E-01  2.1616E-01  4.6060E-02
            -3.1357E-02
 GRADIENT:   5.3347E+02  5.1370E+02  3.6820E+00  1.1884E+02  1.5329E+01  5.7416E+01  6.9021E+00  1.9720E-01  8.3765E+00  9.7977E-01
             8.5057E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1690.62997791032        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     1765             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7561E-01  1.5307E+00  6.5451E-01  6.6673E-01  1.0595E+00  9.4262E-01  8.9837E-01  1.0224E+00  1.1229E+00  9.4708E-01
             8.7687E-01
 PARAMETER:  7.5311E-02  5.2571E-01 -3.2387E-01 -3.0538E-01  1.5779E-01  4.0905E-02 -7.1714E-03  1.2212E-01  2.1592E-01  4.5630E-02
            -3.1397E-02
 GRADIENT:   5.3352E+02  5.1394E+02  3.8629E+00  1.1872E+02  1.5053E+01  5.7421E+01  6.8720E+00  1.5639E-01  8.3015E+00  9.2725E-01
             8.1327E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1690.63011310218        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1946
 NPARAMETR:  9.7562E-01  1.5306E+00  6.5384E-01  6.6678E-01  1.0592E+00  9.4262E-01  8.9849E-01  1.0218E+00  1.1228E+00  9.4681E-01
             8.7688E-01
 PARAMETER:  7.5313E-02  5.2564E-01 -3.2490E-01 -3.0529E-01  1.5748E-01  4.0910E-02 -7.0360E-03  1.2161E-01  2.1587E-01  4.5343E-02
            -3.1381E-02
 GRADIENT:   1.9813E+00 -3.4936E+00 -8.3614E-02 -1.7395E-01  4.6794E-01  1.6916E-01  2.9834E-02  4.1232E-02  8.7983E-02  7.2963E-02
             3.7270E-02

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1690.63035906876        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     2132             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7562E-01  1.5306E+00  6.5352E-01  6.6675E-01  1.0587E+00  9.4263E-01  8.9851E-01  1.0181E+00  1.1225E+00  9.4636E-01
             8.7684E-01
 PARAMETER:  7.5318E-02  5.2564E-01 -3.2538E-01 -3.0534E-01  1.5706E-01  4.0915E-02 -7.0172E-03  1.1794E-01  2.1552E-01  4.4863E-02
            -3.1428E-02
 GRADIENT:   5.3352E+02  5.1393E+02  3.9902E+00  1.1867E+02  1.4804E+01  5.7421E+01  6.8557E+00  1.2434E-01  8.2480E+00  8.9443E-01
             7.9006E-01

0ITERATION NO.:   62    OBJECTIVE VALUE:  -1690.63035906876        NO. OF FUNC. EVALS.:  59
 CUMULATIVE NO. OF FUNC. EVALS.:     2191
 NPARAMETR:  9.7562E-01  1.5306E+00  6.5352E-01  6.6675E-01  1.0587E+00  9.4263E-01  8.9851E-01  1.0181E+00  1.1225E+00  9.4636E-01
             8.7684E-01
 PARAMETER:  7.5318E-02  5.2564E-01 -3.2538E-01 -3.0534E-01  1.5706E-01  4.0915E-02 -7.0172E-03  1.1794E-01  2.1552E-01  4.4863E-02
            -3.1428E-02
 GRADIENT:  -1.8228E-03  1.0050E-01  1.4940E-01 -5.7104E-02  1.5176E-01 -6.0746E-04 -1.7170E-02 -4.1092E-03 -6.5121E-03  1.0302E-02
            -9.4979E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     2191
 NO. OF SIG. DIGITS IN FINAL EST.:  2.3

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.9025E-04 -2.0404E-02 -2.8284E-02  1.7834E-02 -3.6251E-02
 SE:             2.9884E-02  2.3976E-02  1.0780E-02  2.1985E-02  2.1648E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9492E-01  3.9475E-01  8.6952E-03  4.1726E-01  9.4023E-02

 ETASHRINKSD(%)  1.0000E-10  1.9678E+01  6.3886E+01  2.6348E+01  2.7476E+01
 ETASHRINKVR(%)  1.0000E-10  3.5483E+01  8.6958E+01  4.5753E+01  4.7403E+01
 EBVSHRINKSD(%)  3.6259E-01  1.9669E+01  6.7440E+01  2.7966E+01  2.4609E+01
 EBVSHRINKVR(%)  7.2386E-01  3.5470E+01  8.9399E+01  4.8111E+01  4.3162E+01
 RELATIVEINF(%)  9.9110E+01  2.7565E+00  6.8019E-01  2.0026E+00  1.0940E+01
 EPSSHRINKSD(%)  4.6078E+01
 EPSSHRINKVR(%)  7.0924E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1690.6303590687562     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -955.47953250501803     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    29.86
 Elapsed covariance  time in seconds:     5.85
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1690.630       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.76E-01  1.53E+00  6.54E-01  6.67E-01  1.06E+00  9.43E-01  8.99E-01  1.02E+00  1.12E+00  9.46E-01  8.77E-01
 


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
 
         2.77E-02  3.23E-01  2.33E-01  2.13E-01  1.29E-01  6.51E-02  1.35E-01  5.87E-01  2.65E-01  1.41E-01  7.03E-02
 


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
+        7.66E-04
 
 TH 2
+       -9.29E-04  1.04E-01
 
 TH 3
+        6.22E-04 -6.07E-02  5.43E-02
 
 TH 4
+        7.20E-04 -6.79E-02  4.16E-02  4.56E-02
 
 TH 5
+       -2.70E-04  3.18E-02 -9.08E-03 -2.03E-02  1.67E-02
 
 TH 6
+        4.31E-04 -3.65E-03  2.35E-03  2.35E-03 -1.39E-03  4.24E-03
 
 TH 7
+        4.81E-04 -3.37E-02  1.83E-02  2.16E-02 -1.04E-02  4.67E-04  1.83E-02
 
 TH 8
+       -1.01E-03  2.66E-02  2.35E-02 -1.58E-02  2.34E-02 -5.87E-05  5.06E-04  3.44E-01
 
 TH 9
+       -1.21E-03  6.50E-02 -3.16E-02 -4.18E-02  2.35E-02 -3.57E-03 -2.59E-02 -8.28E-03  7.00E-02
 
 TH10
+        2.83E-04  1.48E-02 -2.06E-03 -9.21E-03  1.04E-02 -1.24E-03 -5.48E-03 -1.11E-02  1.36E-02  1.99E-02
 
 TH11
+        1.88E-04 -1.59E-03  2.94E-03  1.11E-03  6.09E-04  6.53E-04  8.05E-04  6.22E-03 -2.37E-04 -2.93E-04  4.94E-03
 
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
+        2.77E-02
 
 TH 2
+       -1.04E-01  3.23E-01
 
 TH 3
+        9.64E-02 -8.07E-01  2.33E-01
 
 TH 4
+        1.22E-01 -9.85E-01  8.37E-01  2.13E-01
 
 TH 5
+       -7.56E-02  7.63E-01 -3.01E-01 -7.34E-01  1.29E-01
 
 TH 6
+        2.39E-01 -1.74E-01  1.55E-01  1.69E-01 -1.65E-01  6.51E-02
 
 TH 7
+        1.28E-01 -7.73E-01  5.81E-01  7.48E-01 -5.94E-01  5.30E-02  1.35E-01
 
 TH 8
+       -6.19E-02  1.40E-01  1.72E-01 -1.26E-01  3.08E-01 -1.54E-03  6.37E-03  5.87E-01
 
 TH 9
+       -1.65E-01  7.61E-01 -5.13E-01 -7.39E-01  6.88E-01 -2.07E-01 -7.24E-01 -5.33E-02  2.65E-01
 
 TH10
+        7.25E-02  3.25E-01 -6.26E-02 -3.06E-01  5.69E-01 -1.35E-01 -2.87E-01 -1.34E-01  3.64E-01  1.41E-01
 
 TH11
+        9.64E-02 -7.00E-02  1.79E-01  7.38E-02  6.71E-02  1.43E-01  8.46E-02  1.51E-01 -1.27E-02 -2.96E-02  7.03E-02
 
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
+        1.50E+03
 
 TH 2
+       -1.34E+02  4.66E+02
 
 TH 3
+       -1.18E+01  1.28E+02  3.42E+02
 
 TH 4
+       -1.43E+02  4.14E+02 -3.42E+02  1.12E+03
 
 TH 5
+        3.18E+01 -2.07E+02 -3.50E+02  3.14E+02  6.68E+02
 
 TH 6
+       -1.45E+02  1.20E+01 -4.26E+01  6.68E+01  3.82E+01  2.82E+02
 
 TH 7
+       -4.98E+01  8.76E+01  4.43E+01  2.40E+00 -5.50E+01  3.88E+01  1.70E+02
 
 TH 8
+        6.95E+00 -1.29E+01 -2.63E+01  2.19E+01  6.46E+00  3.57E+00 -5.51E+00  6.85E+00
 
 TH 9
+        3.31E+01 -3.09E+01 -3.62E+01  2.68E+01 -2.34E+00  2.01E+01  2.17E+01  7.98E+00  5.23E+01
 
 TH10
+       -4.82E+01  7.49E+00 -1.54E+01  1.43E+01 -9.38E+01  1.59E+01  1.19E+00  1.06E+01  8.19E+00  1.01E+02
 
 TH11
+       -4.76E+01  1.10E+01 -1.30E+01  1.64E+01 -2.22E+01 -3.13E+01 -1.48E+01 -1.67E+00 -8.64E+00  1.25E+01  2.23E+02
 
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
 #CPUT: Total CPU Time in Seconds,       35.771
Stop Time:
Wed Sep 29 14:45:38 CDT 2021
