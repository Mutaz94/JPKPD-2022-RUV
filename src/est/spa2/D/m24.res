Thu Sep 30 08:49:15 CDT 2021
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
$DATA ../../../../data/spa2/D/dat24.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m24.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   6482.24854121038        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   2.3558E+02  4.1021E+01  3.6665E+01 -1.9411E+02  2.2592E+02 -1.0843E+03 -3.9267E+02 -1.2160E+02 -7.8908E+02 -4.4788E+02
            -1.4958E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -945.982070535490        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  8.8703E-01  1.1613E+00  9.5477E-01  1.8437E+00  9.1326E-01  3.1978E+00  3.2947E+00  1.0417E+00  3.7843E+00  2.2000E+00
             8.2940E+00
 PARAMETER: -1.9879E-02  2.4957E-01  5.3718E-02  7.1177E-01  9.2692E-03  1.2625E+00  1.2923E+00  1.4086E-01  1.4309E+00  8.8844E-01
             2.2155E+00
 GRADIENT:  -6.7610E+01 -1.5022E+01 -3.5436E+01  5.1439E+01 -2.7548E+01  1.4674E+02  1.6981E+01  4.5537E+00  1.1031E+02  2.9152E+01
             3.9697E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1000.78142187012        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      154
 NPARAMETR:  1.0636E+00  1.5691E+00  6.1244E+00  1.6394E+00  3.7468E+00  3.1758E+00  5.6640E+00  1.1701E+00  7.1241E+00  1.8969E+00
             7.8707E+00
 PARAMETER:  1.6165E-01  5.5049E-01  1.9123E+00  5.9434E-01  1.4209E+00  1.2556E+00  1.8341E+00  2.5713E-01  2.0635E+00  7.4024E-01
             2.1631E+00
 GRADIENT:  -2.2369E+01 -3.5598E-01 -1.5874E+01  3.3600E+01 -4.6504E+00  1.4317E+02  1.0705E+02  1.4965E-01  1.4634E+02  8.3905E+00
             3.7949E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1133.18800280277        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      226
 NPARAMETR:  1.3008E+00  3.6039E-01  8.3621E+00  1.8982E+00  1.2896E+00  2.4567E+00  6.4535E+00  2.2082E+00  1.9543E+00  1.8322E+00
             6.5595E+00
 PARAMETER:  3.6296E-01 -9.2057E-01  2.2237E+00  7.4088E-01  3.5437E-01  9.9882E-01  1.9646E+00  8.9217E-01  7.7001E-01  7.0551E-01
             1.9809E+00
 GRADIENT:   6.9373E+01  1.2419E+01  1.2329E+01  1.8037E+01 -9.5913E+01  5.2179E+01  9.7817E+01  2.5599E+00  5.0585E+00  3.1593E+01
             2.6538E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1186.20937576735        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      297
 NPARAMETR:  1.0694E+00  1.5608E+00  1.2040E+00  7.6139E-01  1.6042E+00  2.0822E+00  2.8864E+00  7.9863E-01  1.2157E+00  6.9893E-01
             4.8925E+00
 PARAMETER:  1.6710E-01  5.4517E-01  2.8569E-01 -1.7262E-01  5.7260E-01  8.3342E-01  1.1600E+00 -1.2486E-01  2.9530E-01 -2.5821E-01
             1.6877E+00
 GRADIENT:  -8.8054E+00 -1.7913E+01 -8.6074E+00  1.9822E+00  1.7559E+01  2.5091E+00 -1.0937E+01  8.9821E-01  1.3401E+01  7.4001E+00
            -3.9271E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1205.20963404344        NO. OF FUNC. EVALS.: 114
 CUMULATIVE NO. OF FUNC. EVALS.:      411
 NPARAMETR:  1.1120E+00  1.7398E+00  1.8625E+00  7.2052E-01  1.9456E+00  2.2295E+00  3.1906E+00  8.9829E-01  9.7099E-01  5.6093E-01
             5.1294E+00
 PARAMETER:  2.0619E-01  6.5375E-01  7.2193E-01 -2.2778E-01  7.6558E-01  9.0180E-01  1.2602E+00 -7.2652E-03  7.0557E-02 -4.7816E-01
             1.7350E+00
 GRADIENT:  -1.8827E+01 -2.2109E+01 -3.8398E+00 -1.6174E+01  2.1742E+01 -5.7625E+01 -3.2838E+01  4.1994E-01  8.7380E+00  3.7525E+00
             6.9005E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1231.05765549783        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      589
 NPARAMETR:  1.1341E+00  1.3749E+00  3.2773E+00  9.9985E-01  1.8742E+00  2.1200E+00  4.6036E+00  1.1250E+00  6.6595E-01  3.0707E-01
             5.2482E+00
 PARAMETER:  2.2587E-01  4.1838E-01  1.2870E+00  9.9848E-02  7.2816E-01  8.5141E-01  1.6268E+00  2.1778E-01 -3.0654E-01 -1.0807E+00
             1.7579E+00
 GRADIENT:   7.4716E+00  4.6525E+00  2.2684E+00 -1.2088E+01 -9.7330E+00 -1.1636E+01  7.9855E+00  5.0265E-01  5.3492E+00  8.5346E-01
            -7.4995E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1233.35624218648        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      769
 NPARAMETR:  1.1243E+00  1.3434E+00  2.6105E+00  9.5085E-01  1.8464E+00  2.1836E+00  4.5013E+00  7.5742E-01  3.3964E-01  1.3784E-01
             5.2673E+00
 PARAMETER:  2.1715E-01  3.9523E-01  1.0595E+00  4.9597E-02  7.1324E-01  8.8099E-01  1.6044E+00 -1.7783E-01 -9.7987E-01 -1.8817E+00
             1.7615E+00
 GRADIENT:   2.9169E+00 -2.0004E+00  5.2265E-01 -6.1772E+00 -1.2638E+00  5.7507E-01  2.6605E+00  4.0691E-01  2.9422E-01  1.6612E-01
            -3.1670E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1233.59135165553        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      944
 NPARAMETR:  1.1167E+00  1.4520E+00  2.2124E+00  9.1235E-01  1.8033E+00  2.1769E+00  4.3059E+00  3.8970E-01  2.1185E-01  6.9917E-02
             5.2862E+00
 PARAMETER:  2.1037E-01  4.7296E-01  8.9406E-01  8.2726E-03  6.8964E-01  8.7790E-01  1.5600E+00 -8.4238E-01 -1.4519E+00 -2.5605E+00
             1.7651E+00
 GRADIENT:  -5.7383E-01  3.0471E-01  1.3746E-01 -8.7784E-02 -7.9261E-01 -3.7182E-01 -7.0997E-01  1.3884E-01 -6.9713E-02  4.4010E-02
            -1.0602E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1233.65774204324        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1124
 NPARAMETR:  1.1181E+00  1.4597E+00  2.1475E+00  9.0528E-01  1.7977E+00  2.1876E+00  4.3049E+00  1.8858E-01  2.1559E-01  3.4345E-02
             5.2884E+00
 PARAMETER:  2.1166E-01  4.7822E-01  8.6431E-01  4.8488E-04  6.8651E-01  8.8282E-01  1.5598E+00 -1.5682E+00 -1.4344E+00 -3.2713E+00
             1.7655E+00
 GRADIENT:   1.1910E-02 -2.4244E-01 -3.1776E-01 -4.7711E+00  6.4496E-01  1.4250E+00  1.5822E+00  3.4511E-02  9.9603E-02  1.0956E-02
             1.6444E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1233.70182494502        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1300
 NPARAMETR:  1.1174E+00  1.4288E+00  2.2059E+00  9.2212E-01  1.7967E+00  2.1767E+00  4.3337E+00  3.6263E-02  2.5752E-01  1.0000E-02
             5.2863E+00
 PARAMETER:  2.1105E-01  4.5686E-01  8.9112E-01  1.8916E-02  6.8593E-01  8.7780E-01  1.5664E+00 -3.2169E+00 -1.2567E+00 -4.8690E+00
             1.7651E+00
 GRADIENT:  -2.4868E-01 -6.1172E-02 -1.0300E-01 -1.2604E-01 -1.4454E-01 -4.2982E-01 -8.6615E-01  1.2335E-03  1.3938E-02  0.0000E+00
            -9.3673E-02

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1233.73464375376        NO. OF FUNC. EVALS.: 166
 CUMULATIVE NO. OF FUNC. EVALS.:     1466             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1187E+00  1.4286E+00  2.2144E+00  9.2361E-01  1.7984E+00  2.2025E+00  4.3576E+00  1.0000E-02  2.5245E-01  1.0000E-02
             5.2854E+00
 PARAMETER:  2.1212E-01  4.5671E-01  8.9500E-01  2.0530E-02  6.8689E-01  8.8960E-01  1.5719E+00 -7.8193E+00 -1.2766E+00 -4.8690E+00
             1.7650E+00
 GRADIENT:   3.2600E+01  1.3916E+01  1.4288E-01  3.1842E+00  2.3962E+00  7.6215E+01  1.0038E+02  0.0000E+00  1.6292E-01  0.0000E+00
             2.0810E+01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1233.73623546119        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     1651             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1187E+00  1.4236E+00  2.2267E+00  9.2641E-01  1.7995E+00  2.2027E+00  4.3659E+00  1.0000E-02  2.5472E-01  1.0000E-02
             5.2852E+00
 PARAMETER:  2.1213E-01  4.5318E-01  9.0050E-01  2.3561E-02  6.8749E-01  8.8969E-01  1.5738E+00 -7.8193E+00 -1.2676E+00 -4.8690E+00
             1.7649E+00
 GRADIENT:   3.2592E+01  1.3815E+01  1.4835E-01  4.2291E+00  2.3697E+00  7.6252E+01  1.0042E+02  0.0000E+00  1.1782E-01  0.0000E+00
             2.0473E+01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1233.73979192072        NO. OF FUNC. EVALS.: 194
 CUMULATIVE NO. OF FUNC. EVALS.:     1845             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1186E+00  1.4193E+00  2.2361E+00  9.2694E-01  1.8005E+00  2.2037E+00  4.3732E+00  1.0000E-02  2.6179E-01  1.0000E-02
             5.2854E+00
 PARAMETER:  2.1206E-01  4.5016E-01  9.0474E-01  2.4133E-02  6.8804E-01  8.9012E-01  1.5755E+00 -7.8193E+00 -1.2402E+00 -4.8690E+00
             1.7649E+00
 GRADIENT:   3.2585E+01  1.3472E+01  1.7225E-01  3.0195E+00  2.3931E+00  7.6460E+01  1.0123E+02  0.0000E+00  1.7835E-01  0.0000E+00
             2.0887E+01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1233.74104436388        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:     2034
 NPARAMETR:  1.1186E+00  1.4178E+00  2.2431E+00  9.2823E-01  1.8011E+00  2.2027E+00  4.3783E+00  1.0000E-02  2.6511E-01  1.0000E-02
             5.2854E+00
 PARAMETER:  2.1211E-01  4.4908E-01  9.0785E-01  2.5522E-02  6.8841E-01  8.8967E-01  1.5767E+00 -7.8193E+00 -1.2276E+00 -4.8690E+00
             1.7650E+00
 GRADIENT:   2.3737E-01  1.5075E-01 -8.6914E-02  2.2793E-01 -1.5644E-01  4.0909E+00  3.4109E-01  0.0000E+00 -3.3311E-02  0.0000E+00
            -3.3514E-01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1233.74258163223        NO. OF FUNC. EVALS.: 195
 CUMULATIVE NO. OF FUNC. EVALS.:     2229             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1187E+00  1.4142E+00  2.2518E+00  9.2916E-01  1.8021E+00  2.2026E+00  4.3854E+00  1.0000E-02  2.7178E-01  1.0000E-02
             5.2859E+00
 PARAMETER:  2.1212E-01  4.4654E-01  9.1173E-01  2.6520E-02  6.8895E-01  8.8965E-01  1.5783E+00 -7.8193E+00 -1.2028E+00 -4.8690E+00
             1.7650E+00
 GRADIENT:   3.2635E+01  1.3256E+01  1.7822E-01  2.4397E+00  2.4420E+00  7.6205E+01  1.0201E+02  0.0000E+00  2.3072E-01  0.0000E+00
             2.1212E+01

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1233.74309863638        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     2414
 NPARAMETR:  1.1186E+00  1.4114E+00  2.2592E+00  9.2989E-01  1.8032E+00  2.2026E+00  4.3909E+00  1.0000E-02  2.7581E-01  1.0000E-02
             5.2861E+00
 PARAMETER:  2.1211E-01  4.4462E-01  9.1502E-01  2.7315E-02  6.8958E-01  8.8965E-01  1.5795E+00 -7.8193E+00 -1.1880E+00 -4.8690E+00
             1.7651E+00
 GRADIENT:   2.5653E-01 -9.1190E-02 -9.4635E-02 -1.0163E+00 -7.0453E-03  4.0611E+00  8.9973E-01  0.0000E+00  3.6737E-02  0.0000E+00
             2.2785E-01

0ITERATION NO.:   85    OBJECTIVE VALUE:  -1233.74420696086        NO. OF FUNC. EVALS.: 195
 CUMULATIVE NO. OF FUNC. EVALS.:     2609             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1186E+00  1.4102E+00  2.2668E+00  9.3160E-01  1.8039E+00  2.2027E+00  4.3926E+00  1.0000E-02  2.7502E-01  1.0000E-02
             5.2856E+00
 PARAMETER:  2.1210E-01  4.4375E-01  9.1838E-01  2.9143E-02  6.8992E-01  8.8968E-01  1.5799E+00 -7.8193E+00 -1.1909E+00 -4.8690E+00
             1.7650E+00
 GRADIENT:   3.2620E+01  1.3201E+01  2.0638E-01  3.1904E+00  2.3880E+00  7.6215E+01  1.0207E+02  0.0000E+00  2.0077E-01  0.0000E+00
             2.0932E+01

0ITERATION NO.:   90    OBJECTIVE VALUE:  -1233.74459248264        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:     2798
 NPARAMETR:  1.1186E+00  1.4089E+00  2.2719E+00  9.3230E-01  1.8045E+00  2.2027E+00  4.3951E+00  1.0000E-02  2.7653E-01  1.0000E-02
             5.2855E+00
 PARAMETER:  2.1209E-01  4.4282E-01  9.2061E-01  2.9904E-02  6.9028E-01  8.8968E-01  1.5805E+00 -7.8193E+00 -1.1854E+00 -4.8690E+00
             1.7650E+00
 GRADIENT:   2.4031E-01  8.2037E-02 -5.1756E-02  2.8679E-01 -1.4565E-01  4.0892E+00  4.6029E-01  0.0000E+00 -2.7029E-02  0.0000E+00
            -2.8985E-01

0ITERATION NO.:   93    OBJECTIVE VALUE:  -1233.74479365720        NO. OF FUNC. EVALS.:  97
 CUMULATIVE NO. OF FUNC. EVALS.:     2895
 NPARAMETR:  1.1186E+00  1.4077E+00  2.2737E+00  9.3206E-01  1.8049E+00  2.2027E+00  4.3979E+00  1.0000E-02  2.7982E-01  1.0000E-02
             5.2860E+00
 PARAMETER:  2.1210E-01  4.4198E-01  9.2142E-01  2.9641E-02  6.9049E-01  8.8967E-01  1.5811E+00 -7.8193E+00 -1.1736E+00 -4.8690E+00
             1.7651E+00
 GRADIENT:   2.5331E-01 -5.0388E-02 -6.3902E-02 -5.6377E-01 -5.8386E-02  4.0699E+00  7.9204E-01  0.0000E+00  1.9060E-02  0.0000E+00
             7.6584E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     2895
 NO. OF SIG. DIGITS IN FINAL EST.:  2.0
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         6.6722E-03  1.0521E-02 -6.2254E-06 -2.5007E-02  2.2063E-05
 SE:             2.9280E-02  2.7389E-02  3.5133E-05  5.6110E-03  1.2367E-04
 N:                     100         100         100         100         100

 P VAL.:         8.1974E-01  7.0088E-01  8.5935E-01  8.3249E-06  8.5841E-01

 ETASHRINKSD(%)  1.9077E+00  8.2418E+00  9.9882E+01  8.1202E+01  9.9586E+01
 ETASHRINKVR(%)  3.7790E+00  1.5804E+01  1.0000E+02  9.6467E+01  9.9998E+01
 EBVSHRINKSD(%)  2.0247E+00  4.7376E+00  9.9867E+01  8.6841E+01  9.9551E+01
 EBVSHRINKVR(%)  4.0083E+00  9.2507E+00  1.0000E+02  9.8268E+01  9.9998E+01
 RELATIVEINF(%)  9.5818E+01  4.9216E+01  2.9591E-05  8.8317E-01  3.5141E-04
 EPSSHRINKSD(%)  1.3133E+01
 EPSSHRINKVR(%)  2.4542E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          600
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1102.7262398456071     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1233.7447936571975     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -131.01855381159044     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    63.88
0R MATRIX ALGORITHMICALLY SINGULAR
0COVARIANCE MATRIX UNOBTAINABLE
0R MATRIX IS OUTPUT
0S MATRIX ALGORITHMICALLY SINGULAR
0S MATRIX IS OUTPUT
0T MATRIX - EQUAL TO RS*RMAT, WHERE S* IS A PSEUDO INVERSE OF S - IS OUTPUT
 Elapsed covariance  time in seconds:    10.28
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1233.745       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.12E+00  1.41E+00  2.27E+00  9.32E-01  1.80E+00  2.20E+00  4.40E+00  1.00E-02  2.80E-01  1.00E-02  5.29E+00
 


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
+        1.23E+02
 
 TH 2
+       -1.20E+00  5.09E+00
 
 TH 3
+        3.55E+00 -8.30E-01  2.26E-01
 
 TH 4
+       -5.69E+01  4.69E+01 -8.89E+00  4.49E+02
 
 TH 5
+       -1.30E+01  1.53E+00 -5.93E-01  1.88E+01  1.75E+00
 
 TH 6
+        1.04E+01 -2.34E-01  3.18E-01 -5.99E+00 -1.12E+00  8.72E-01
 
 TH 7
+        4.62E+00 -3.30E+00  6.42E-01 -3.18E+01 -1.39E+00  4.73E-01  2.26E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        7.24E+00 -6.17E+00  1.16E+00 -5.90E+01 -2.45E+00  7.67E-01  4.18E+00  0.00E+00  7.75E+00
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -7.31E-01 -2.24E+00  3.26E-01 -2.01E+01 -5.33E-01  2.37E-03  1.41E+00  0.00E+00  2.66E+00  0.00E+00  1.66E+00
 
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
+        1.74E+02
 
 TH 2
+       -1.26E+00  2.63E+01
 
 TH 3
+        7.09E-01  1.26E+00  4.27E+00
 
 TH 4
+       -7.16E+00  4.72E+01 -8.18E+00  4.39E+02
 
 TH 5
+       -2.63E+00 -7.49E+00 -1.47E+01  1.56E+01  6.35E+01
 
 TH 6
+       -2.85E-01 -7.73E-02  1.53E-01  2.03E+00 -1.05E+00  3.82E+01
 
 TH 7
+        3.98E-01  1.77E+00 -5.91E-01 -3.08E+01  2.07E+00 -4.05E-01  7.67E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        4.71E-01 -2.27E+00 -1.09E+00 -5.68E+01  4.60E+00 -4.22E-01  3.19E+00  0.00E+00  1.72E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -4.84E+00 -2.76E+00  1.67E-01 -2.04E+01  4.82E-02  6.04E-01  1.23E+00  0.00E+00  3.10E+00  0.00E+00  2.66E+01
 
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
+        1.76E+02
 
 TH 2
+        4.97E+01  2.64E+01
 
 TH 3
+        3.11E+00  1.38E+00  2.83E+00
 
 TH 4
+        7.15E+01  5.03E+01 -6.21E-02  4.45E+02
 
 TH 5
+       -2.24E+01 -4.75E+00 -1.12E+01 -1.18E+01  5.41E+01
 
 TH 6
+        4.36E+01  9.29E+00  4.30E-02 -1.48E+01 -8.48E+00  6.28E+01
 
 TH 7
+        1.40E+00  4.71E+00 -7.68E-01 -3.38E+01  8.76E+00  8.61E+00  1.52E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -1.11E+00 -2.72E+00 -1.08E+00 -5.40E+01  6.23E+00 -1.86E-01  3.44E+00  0.00E+00  1.14E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -1.14E+01 -2.23E+00  2.56E-01 -5.85E+00  4.76E+00 -2.18E+01  5.35E+00  0.00E+00  2.37E-01  0.00E+00  9.71E+02
 
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
 #CPUT: Total CPU Time in Seconds,       74.226
Stop Time:
Thu Sep 30 08:50:31 CDT 2021
