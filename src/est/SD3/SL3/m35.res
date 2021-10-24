Sun Oct 24 00:12:22 CDT 2021
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
$DATA ../../../../data/SD3/SL3/dat35.csv ignore=@
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

NM-TRAN MESSAGES
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.

License Registered to: University of Minnesota
Expiration Date:    14 APR 2022
Current Date:       24 OCT 2021
Days until program expires : 175
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

 #PARA: PARAFILE=../../../../rmpi.pnm, PROTOCOL=MPI, NODES= 10

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
 RAW OUTPUT FILE (FILE): m35.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2034.39529316478        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.4631E+02 -9.6966E+00 -7.8884E+01  9.7743E+01  1.4287E+02  6.1846E+01  3.5048E+00  1.4614E+01  8.6558E+00 -2.1610E+01
            -2.1922E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2047.46347137338        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  1.1341E+00  1.0088E+00  1.1309E+00  9.7826E-01  9.1822E-01  9.4642E-01  9.6418E-01  8.8733E-01  9.6672E-01  9.6988E-01
             1.5271E+00
 PARAMETER:  2.2582E-01  1.0873E-01  2.2301E-01  7.8022E-02  1.4681E-02  4.4931E-02  6.3526E-02 -1.9540E-02  6.6154E-02  6.9419E-02
             5.2340E-01
 GRADIENT:   6.2586E+02  2.4059E+00  1.7599E+01 -3.2765E+01 -5.4195E+01 -2.0649E+00  3.2330E+00  7.8545E+00 -3.2961E+00  2.3515E+00
             1.4209E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2060.17883226496        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  1.1022E+00  1.0654E+00  1.1502E+00  9.4727E-01  9.4354E-01  8.7413E-01  8.1194E-01  2.7609E-01  1.1299E+00  1.1277E+00
             1.4129E+00
 PARAMETER:  1.9732E-01  1.6336E-01  2.3990E-01  4.5825E-02  4.1885E-02 -3.4530E-02 -1.0832E-01 -1.1870E+00  2.2211E-01  2.2022E-01
             4.4562E-01
 GRADIENT:   5.8389E+02  3.6805E+01  2.1412E+01 -4.8028E+00 -6.5189E+01 -2.7610E+01  6.3857E+00  5.8790E-01  2.1135E+01  1.0199E+01
             1.0110E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2066.88353283582        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      229
 NPARAMETR:  1.0192E+00  1.0279E+00  1.2963E+00  9.6701E-01  9.9795E-01  8.5929E-01  8.6246E-01  2.5441E-01  1.0706E+00  1.1995E+00
             1.2659E+00
 PARAMETER:  1.1900E-01  1.2750E-01  3.5950E-01  6.6453E-02  9.7947E-02 -5.1651E-02 -4.7967E-02 -1.2688E+00  1.6823E-01  2.8192E-01
             3.3576E-01
 GRADIENT:   2.4040E+02  2.4497E+01  2.1313E+01 -6.8592E+00 -4.5015E+01 -2.5180E+01  6.7698E+00  3.2966E-01  1.0547E+01  8.5386E+00
             2.2676E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2075.96667962311        NO. OF FUNC. EVALS.: 157
 CUMULATIVE NO. OF FUNC. EVALS.:      386
 NPARAMETR:  1.0562E+00  9.4648E-01  1.3739E+00  1.0600E+00  1.0335E+00  9.9394E-01  7.9216E-01  2.7699E-01  1.0330E+00  1.1981E+00
             1.2368E+00
 PARAMETER:  1.5466E-01  4.4995E-02  4.1766E-01  1.5828E-01  1.3299E-01  9.3925E-02 -1.3299E-01 -1.1838E+00  1.3243E-01  2.8072E-01
             3.1255E-01
 GRADIENT:  -4.2099E+00  8.5029E+00 -5.8587E-01  5.5894E+00  6.9920E-01  6.1486E+00  8.2538E-01  2.9290E-01  1.9042E+00  6.0700E-01
            -3.2143E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2077.72561818647        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      561
 NPARAMETR:  1.0486E+00  5.7689E-01  1.5367E+00  1.3010E+00  9.5317E-01  9.7652E-01  7.2979E-01  1.2504E-01  8.9629E-01  1.1808E+00
             1.2386E+00
 PARAMETER:  1.4742E-01 -4.5011E-01  5.2963E-01  3.6313E-01  5.2037E-02  7.6239E-02 -2.1500E-01 -1.9791E+00 -9.4874E-03  2.6615E-01
             3.1402E-01
 GRADIENT:  -1.2501E+01  1.0164E+01  4.2697E+00  2.2292E+01 -8.7296E+00  7.5286E-01 -6.3384E-02  2.6268E-02  8.8283E-01 -1.6207E+00
            -8.6972E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2078.31827395227        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      738
 NPARAMETR:  1.0531E+00  3.5059E-01  1.6446E+00  1.4364E+00  9.2636E-01  9.7663E-01  6.5799E-01  3.0433E-02  8.2067E-01  1.1898E+00
             1.2400E+00
 PARAMETER:  1.5169E-01 -9.4814E-01  5.9749E-01  4.6211E-01  2.3508E-02  7.6351E-02 -3.1857E-01 -3.3922E+00 -9.7636E-02  2.7379E-01
             3.1507E-01
 GRADIENT:   4.5705E+00  2.0347E+00  1.0699E+00  4.9294E+00 -8.9292E-01  1.9371E+00  8.1759E-02  5.9793E-04  6.5706E-01 -1.4937E+00
            -1.1426E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2078.37712754005        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      920
 NPARAMETR:  1.0507E+00  3.3886E-01  1.6432E+00  1.4394E+00  9.2400E-01  9.7148E-01  4.1054E-01  1.0000E-02  8.1999E-01  1.1987E+00
             1.2407E+00
 PARAMETER:  1.4944E-01 -9.8217E-01  5.9667E-01  4.6423E-01  2.0953E-02  7.1069E-02 -7.9027E-01 -4.5707E+00 -9.8469E-02  2.8124E-01
             3.1568E-01
 GRADIENT:  -2.6277E-01  8.7985E-01  3.0765E-01 -2.3327E+00 -2.0516E-01 -1.1901E-02  2.5953E-02  0.0000E+00  1.6394E-01  9.3622E-02
            -5.4500E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2078.38257939907        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1096
 NPARAMETR:  1.0493E+00  3.0980E-01  1.6478E+00  1.4591E+00  9.1759E-01  9.7080E-01  1.6440E-01  1.0000E-02  8.1300E-01  1.1980E+00
             1.2409E+00
 PARAMETER:  1.4815E-01 -1.0718E+00  5.9944E-01  4.7780E-01  1.3998E-02  7.0363E-02 -1.7054E+00 -5.8084E+00 -1.0702E-01  2.8064E-01
             3.1583E-01
 GRADIENT:  -2.3599E+00  1.1909E+00 -3.7434E-01  1.4404E+00 -2.1619E-01 -1.5836E-01  5.8483E-03  0.0000E+00  7.0415E-01  2.2153E-01
            -1.2418E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2078.39916989065        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     1282             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0511E+00  3.0164E-01  1.6505E+00  1.4600E+00  9.1691E-01  9.7114E-01  9.0774E-02  1.0000E-02  8.1006E-01  1.1973E+00
             1.2411E+00
 PARAMETER:  1.4981E-01 -1.0985E+00  6.0109E-01  4.7846E-01  1.3256E-02  7.0711E-02 -2.2994E+00 -6.1030E+00 -1.1064E-01  2.8010E-01
             3.1596E-01
 GRADIENT:   4.5988E+02  2.9827E+01  4.9051E+00  5.6381E+02  6.0549E+00  3.3431E+01  4.9448E-02  0.0000E+00  7.4882E+00  2.4219E+00
             3.1671E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -2078.40047100386        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1465
 NPARAMETR:  1.0512E+00  3.0194E-01  1.6521E+00  1.4606E+00  9.1622E-01  9.7123E-01  3.0190E-02  1.0000E-02  8.0830E-01  1.1971E+00
             1.2411E+00
 PARAMETER:  1.4994E-01 -1.0975E+00  6.0203E-01  4.7888E-01  1.2507E-02  7.0811E-02 -3.4003E+00 -6.1030E+00 -1.1283E-01  2.7988E-01
             3.1600E-01
 GRADIENT:   2.1850E+00  3.1330E-01  3.1399E-01 -6.2908E+00 -5.4945E-01  7.8671E-02  3.8459E-04  0.0000E+00 -2.6728E-01  5.9075E-02
            -9.9446E-02

0ITERATION NO.:   55    OBJECTIVE VALUE:  -2078.40088881792        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     1653
 NPARAMETR:  1.0513E+00  3.0147E-01  1.6512E+00  1.4607E+00  9.1629E-01  9.7125E-01  1.1089E-02  1.0000E-02  8.0883E-01  1.1969E+00
             1.2412E+00
 PARAMETER:  1.5005E-01 -1.0991E+00  6.0152E-01  4.7890E-01  1.2573E-02  7.0834E-02 -4.4018E+00 -6.1030E+00 -1.1216E-01  2.7973E-01
             3.1607E-01
 GRADIENT:   2.4755E+00  1.5492E-01  4.8367E-02 -6.7796E+00  4.3566E-03  9.0308E-02  7.7142E-05  0.0000E+00  4.9805E-02  2.6150E-02
            -5.5921E-03

0ITERATION NO.:   60    OBJECTIVE VALUE:  -2078.40092333026        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1833
 NPARAMETR:  1.0514E+00  3.0000E-01  1.6442E+00  1.4590E+00  9.1681E-01  9.7131E-01  1.0000E-02  1.0000E-02  8.0964E-01  1.1954E+00
             1.2414E+00
 PARAMETER:  1.5010E-01 -1.0987E+00  6.0127E-01  4.7884E-01  1.2422E-02  7.0842E-02 -4.5392E+00 -6.1030E+00 -1.1228E-01  2.7963E-01
             3.1606E-01
 GRADIENT:  -1.0745E-02  2.4759E-02  9.1212E-02  2.3773E-01 -5.6909E-02 -1.3624E-03  0.0000E+00  0.0000E+00 -2.7141E-02  1.6958E-02
            -1.0013E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1833
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         4.6937E-04 -1.2099E-04 -1.2083E-04 -4.3697E-03 -2.2969E-02
 SE:             2.9835E-02  5.6278E-05  1.1666E-04  2.9315E-02  2.5153E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8745E-01  3.1571E-02  3.0031E-01  8.8151E-01  3.6114E-01

 ETASHRINKSD(%)  4.7653E-02  9.9811E+01  9.9609E+01  1.7910E+00  1.5734E+01
 ETASHRINKVR(%)  9.5283E-02  1.0000E+02  9.9998E+01  3.5498E+00  2.8993E+01
 EBVSHRINKSD(%)  4.8749E-01  9.9825E+01  9.9608E+01  2.0604E+00  1.2788E+01
 EBVSHRINKVR(%)  9.7260E-01  1.0000E+02  9.9998E+01  4.0783E+00  2.3940E+01
 RELATIVEINF(%)  9.7982E+01  2.0261E-05  2.5821E-04  6.8719E+00  1.3181E+01
 EPSSHRINKSD(%)  3.1158E+01
 EPSSHRINKVR(%)  5.2608E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2078.4009233302613     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1159.4623901255886     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     8.99
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2078.401       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.05E+00  3.02E-01  1.65E+00  1.46E+00  9.16E-01  9.71E-01  1.00E-02  1.00E-02  8.09E-01  1.20E+00  1.24E+00
 


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
 
 Elapsed finaloutput time in seconds:     0.00
1THERE ARE ERROR MESSAGES IN FILE PRDERR                                                                  
 #CPUT: Total CPU Time in Seconds,       58.295
Stop Time:
Sun Oct 24 00:12:34 CDT 2021
