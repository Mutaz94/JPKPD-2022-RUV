Sat Sep 18 10:12:58 CDT 2021
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
$DATA ../../../../data/spa/A3/dat2.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m2.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   224.216985294171        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.2421E+02  8.2146E+01  9.0467E+01 -1.6031E+01  1.5101E+02  8.5616E+00 -7.5239E+01 -2.0043E+01 -1.3431E+02 -1.8475E+02
            -3.3147E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1186.67020663668        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0742E+00  9.8091E-01  9.3429E-01  1.1688E+00  9.3256E-01  8.9499E-01  1.0063E+00  9.7704E-01  1.0740E+00  1.0781E+00
             5.4343E+00
 PARAMETER:  1.7157E-01  8.0724E-02  3.2028E-02  2.5597E-01  3.0174E-02 -1.0944E-02  1.0626E-01  7.6772E-02  1.7141E-01  1.7516E-01
             1.7927E+00
 GRADIENT:   8.6030E+01 -2.1459E+00 -1.0261E+01  9.9544E+00 -1.6632E+01  1.7610E-01  1.3764E+01  6.7031E+00  3.0468E+01  2.5837E+01
             1.7735E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1217.38717527678        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  1.0224E+00  6.9894E-01  4.6157E-01  1.2367E+00  5.1317E-01  9.1272E-01  6.0178E-01  7.3581E-02  1.2123E+00  4.6582E-01
             4.7563E+00
 PARAMETER:  1.2215E-01 -2.5819E-01 -6.7312E-01  3.1241E-01 -5.6715E-01  8.6702E-03 -4.0787E-01 -2.5094E+00  2.9249E-01 -6.6396E-01
             1.6595E+00
 GRADIENT:  -4.5374E+00  1.2395E+01 -1.3923E+01  5.0715E+01  3.0400E-01 -1.0306E+00  3.8646E+00  1.1678E-01  3.3507E+01  1.1424E+01
             1.3407E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1242.33884517171        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      230
 NPARAMETR:  1.0053E+00  4.6582E-01  2.6974E-01  1.2081E+00  3.1801E-01  9.6612E-01  4.3672E-01  2.2603E-02  1.0840E+00  1.8730E-01
             3.6350E+00
 PARAMETER:  1.0526E-01 -6.6396E-01 -1.2103E+00  2.8902E-01 -1.0457E+00  6.5529E-02 -7.2845E-01 -3.6897E+00  1.8068E-01 -1.5750E+00
             1.3906E+00
 GRADIENT:   1.4720E+00  5.1903E+01  2.3032E+01  9.2548E+01 -6.3204E+01  1.6846E+00 -2.7949E+00 -5.6204E-03 -1.9414E+01 -2.5543E+00
            -3.9159E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1249.44329961311        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      300
 NPARAMETR:  9.9799E-01  3.0681E-01  2.5627E-01  1.1844E+00  2.7717E-01  9.7339E-01  5.9912E-01  1.0000E-02  1.1530E+00  1.2574E-01
             3.7090E+00
 PARAMETER:  9.7987E-02 -1.0815E+00 -1.2615E+00  2.6920E-01 -1.1831E+00  7.3033E-02 -4.1229E-01 -4.7693E+00  2.4233E-01 -1.9735E+00
             1.4108E+00
 GRADIENT:  -8.1586E+00  1.2225E+01  2.2797E+01  1.7639E+01 -4.0508E+01  5.8908E+00 -1.2735E+00  0.0000E+00  1.1020E+00 -1.0610E+00
            -7.1983E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1250.84331936245        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      371
 NPARAMETR:  9.9755E-01  2.1412E-01  2.6513E-01  1.2080E+00  2.7021E-01  9.4591E-01  1.4661E+00  1.0000E-02  1.1305E+00  1.0873E-01
             3.7119E+00
 PARAMETER:  9.7551E-02 -1.4412E+00 -1.2275E+00  2.8901E-01 -1.2085E+00  4.4389E-02  4.8261E-01 -6.0138E+00  2.2262E-01 -2.1189E+00
             1.4115E+00
 GRADIENT:  -1.3520E-01  3.3810E+00  8.4558E+00  9.0059E-01 -1.4455E+01 -9.9780E-01  2.3644E-01  0.0000E+00  4.4360E+00 -5.0789E-01
            -3.2280E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1252.19260514638        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      442
 NPARAMETR:  9.8710E-01  8.1198E-02  2.8118E-01  1.2678E+00  2.6926E-01  9.4106E-01  4.0476E+00  1.0000E-02  1.0466E+00  1.1903E-01
             3.7241E+00
 PARAMETER:  8.7012E-02 -2.4109E+00 -1.1688E+00  3.3730E-01 -1.2121E+00  3.9256E-02  1.4981E+00 -8.2513E+00  1.4557E-01 -2.0283E+00
             1.4148E+00
 GRADIENT:  -2.1467E+00  1.1212E+00  4.8393E+00  1.8549E+00 -8.8569E+00 -8.4811E-02  5.8378E-01  0.0000E+00 -1.9631E+00 -3.6728E-01
            -1.6573E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1252.32812564273        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      515
 NPARAMETR:  9.8470E-01  4.4465E-02  2.8053E-01  1.2778E+00  2.6669E-01  9.3913E-01  6.1356E+00  1.0000E-02  1.0515E+00  1.5161E-01
             3.7155E+00
 PARAMETER:  8.4581E-02 -3.0131E+00 -1.1711E+00  3.4510E-01 -1.2217E+00  3.7194E-02  1.9141E+00 -9.5489E+00  1.5019E-01 -1.7864E+00
             1.4125E+00
 GRADIENT:  -9.1932E-02  5.7638E-01 -1.4218E-01  3.6381E-01  6.2689E-02 -5.3491E-01  9.8420E-01  0.0000E+00  3.8552E-01 -4.4720E-01
            -5.6676E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1252.42398117091        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      590
 NPARAMETR:  9.8350E-01  2.8134E-02  2.7977E-01  1.2812E+00  2.6462E-01  9.4001E-01  7.8657E+00  1.0000E-02  1.0536E+00  2.0444E-01
             3.6971E+00
 PARAMETER:  8.3361E-02 -3.4708E+00 -1.1738E+00  3.4778E-01 -1.2295E+00  3.8139E-02  2.1625E+00 -1.0456E+01  1.5222E-01 -1.4875E+00
             1.4076E+00
 GRADIENT:   1.9851E-01  1.3229E+00  1.2472E-01 -1.3882E+00  2.7600E-01 -6.5722E-01  2.6100E+00  0.0000E+00  6.1054E-01 -5.2207E-01
            -5.5595E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1252.49494178441        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      661
 NPARAMETR:  9.8286E-01  1.8462E-02  2.7830E-01  1.2818E+00  2.6229E-01  9.4218E-01  9.7040E+00  1.0000E-02  1.0563E+00  2.7435E-01
             3.6741E+00
 PARAMETER:  8.2712E-02 -3.8920E+00 -1.1790E+00  3.4824E-01 -1.2383E+00  4.0443E-02  2.3725E+00 -1.1280E+01  1.5477E-01 -1.1933E+00
             1.4013E+00
 GRADIENT:  -7.0430E-01  1.0351E+00 -2.2992E-01 -1.9151E+00  1.2241E+00  3.4285E-02  2.0807E+00  0.0000E+00  1.0817E+00 -5.8509E-02
             2.3813E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1252.50474787005        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      732
 NPARAMETR:  9.8326E-01  1.6518E-02  2.7856E-01  1.2829E+00  2.6199E-01  9.4106E-01  1.0191E+01  1.0000E-02  1.0483E+00  2.9391E-01
             3.6474E+00
 PARAMETER:  8.3118E-02 -4.0033E+00 -1.1781E+00  3.4911E-01 -1.2394E+00  3.9249E-02  2.4215E+00 -1.1502E+01  1.4718E-01 -1.1245E+00
             1.3940E+00
 GRADIENT:   4.9295E-01  5.0489E-01  2.6258E-01 -4.8385E-01 -4.1037E-01 -3.1681E-01  1.0073E+00  0.0000E+00 -5.8374E-01 -5.0486E-04
            -1.2069E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1252.51447985832        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      805
 NPARAMETR:  9.8325E-01  1.7159E-02  2.7874E-01  1.2829E+00  2.6219E-01  9.4107E-01  9.9251E+00  1.0000E-02  1.0484E+00  2.9198E-01
             3.6498E+00
 PARAMETER:  8.3107E-02 -3.9652E+00 -1.1775E+00  3.4915E-01 -1.2387E+00  3.9265E-02  2.3951E+00 -1.1414E+01  1.4729E-01 -1.1311E+00
             1.3947E+00
 GRADIENT:   1.9529E-01  1.6216E-01  1.1956E-01 -1.1106E-01 -1.4448E-01 -1.5558E-01  3.3599E-01  0.0000E+00 -3.3166E-01 -1.4635E-02
            -5.6801E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1252.66572908265        NO. OF FUNC. EVALS.: 138
 CUMULATIVE NO. OF FUNC. EVALS.:      943
 NPARAMETR:  9.8506E-01  2.6926E-02  2.9360E-01  1.3021E+00  2.7252E-01  9.3967E-01  7.3371E+00  1.0000E-02  1.0391E+00  2.8034E-01
             3.6631E+00
 PARAMETER:  8.4943E-02 -3.5146E+00 -1.1255E+00  3.6399E-01 -1.2001E+00  3.7770E-02  2.0929E+00 -1.0366E+01  1.3840E-01 -1.1718E+00
             1.3983E+00
 GRADIENT:   7.0087E-01 -3.6329E-02  3.0890E-01  5.5532E-01 -4.6042E-01 -6.3671E-02 -3.1077E-02  0.0000E+00  1.0058E+00 -6.3073E-02
            -1.0771E+00

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1252.67040586273        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1118
 NPARAMETR:  9.8525E-01  3.1479E-02  2.9419E-01  1.3010E+00  2.7329E-01  9.3974E-01  6.5532E+00  1.0000E-02  1.0347E+00  2.7903E-01
             3.6693E+00
 PARAMETER:  8.5142E-02 -3.3584E+00 -1.1235E+00  3.6314E-01 -1.1972E+00  3.7850E-02  1.9799E+00 -9.9910E+00  1.3414E-01 -1.1764E+00
             1.4000E+00
 GRADIENT:   1.1021E-01 -5.4683E-02 -7.2395E-02  8.4164E-03  2.2477E-01 -1.6662E-02 -7.1889E-02  0.0000E+00  1.0200E-01 -2.0390E-02
            -4.6510E-02

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1252.67656668758        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1297
 NPARAMETR:  9.8637E-01  4.4977E-02  2.9270E-01  1.2949E+00  2.7333E-01  9.4014E-01  5.7380E+00  1.0000E-02  1.0374E+00  2.7730E-01
             3.6700E+00
 PARAMETER:  8.6280E-02 -3.0016E+00 -1.1286E+00  3.5840E-01 -1.1971E+00  3.8279E-02  1.8471E+00 -9.2460E+00  1.3670E-01 -1.1827E+00
             1.4002E+00
 GRADIENT:   8.2974E-02 -1.4433E-02 -3.2117E-01 -3.5471E-01  5.5589E-01  3.7920E-02 -9.0381E-03  0.0000E+00  4.7512E-03 -5.1404E-03
             1.1633E-01

0ITERATION NO.:   73    OBJECTIVE VALUE:  -1252.67672161174        NO. OF FUNC. EVALS.:  92
 CUMULATIVE NO. OF FUNC. EVALS.:     1389
 NPARAMETR:  9.8630E-01  4.4204E-02  2.9292E-01  1.2956E+00  2.7336E-01  9.4003E-01  5.8026E+00  1.0000E-02  1.0372E+00  2.7876E-01
             3.6689E+00
 PARAMETER:  8.6207E-02 -3.0189E+00 -1.1279E+00  3.5897E-01 -1.1970E+00  3.8161E-02  1.8583E+00 -9.2834E+00  1.3648E-01 -1.1774E+00
             1.3999E+00
 GRADIENT:   1.1217E-02 -1.2605E-03 -1.3763E-02 -1.2235E-02  3.5110E-02  5.3462E-04 -1.0198E-03  0.0000E+00  5.2176E-03  1.9918E-03
             1.9038E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1389
 NO. OF SIG. DIGITS IN FINAL EST.:  3.2
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -8.8370E-04  1.8779E-04  1.5282E-04 -1.5328E-02  2.5214E-03
 SE:             2.8320E-02  3.4975E-03  2.3575E-04  2.5827E-02  9.6257E-03
 N:                     100         100         100         100         100

 P VAL.:         9.7511E-01  9.5718E-01  5.1684E-01  5.5284E-01  7.9336E-01

 ETASHRINKSD(%)  5.1229E+00  8.8283E+01  9.9210E+01  1.3477E+01  6.7753E+01
 ETASHRINKVR(%)  9.9834E+00  9.8627E+01  9.9994E+01  2.5137E+01  8.9601E+01
 EBVSHRINKSD(%)  4.7447E+00  8.8783E+01  9.9182E+01  1.2574E+01  6.7581E+01
 EBVSHRINKVR(%)  9.2642E+00  9.8742E+01  9.9993E+01  2.3567E+01  8.9490E+01
 RELATIVEINF(%)  6.7937E+01  1.2795E-01  2.3897E-04  1.4077E+01  3.1003E-01
 EPSSHRINKSD(%)  2.3192E+01
 EPSSHRINKVR(%)  4.1005E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1252.6767216117421     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -517.52589504800392     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    14.58
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.87
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1252.677       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.86E-01  4.42E-02  2.93E-01  1.30E+00  2.73E-01  9.40E-01  5.80E+00  1.00E-02  1.04E+00  2.79E-01  3.67E+00
 


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
+        1.18E+03
 
 TH 2
+       -1.49E+02  4.44E+02
 
 TH 3
+       -1.47E+02  7.68E+02  8.18E+03
 
 TH 4
+       -5.80E+01  1.73E+02 -2.64E+02  4.69E+02
 
 TH 5
+        4.14E+02 -1.52E+03 -1.18E+04 -5.01E+02  2.02E+04
 
 TH 6
+        7.06E+00  2.64E+00  4.14E+01 -1.60E+01  1.08E+01  1.85E+02
 
 TH 7
+       -2.84E-02  3.42E+00 -1.31E-01 -1.36E-01  1.05E+00  2.80E-02  5.13E-02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.12E+01 -1.96E+01  5.07E+01 -1.16E+01  1.16E+02  5.50E+00  1.40E-01  0.00E+00  1.04E+02
 
 TH10
+       -1.19E+01 -1.20E+01 -1.52E+02 -1.07E+01  2.68E+02 -7.85E-01  1.12E-01  0.00E+00  2.80E+00  3.42E+01
 
 TH11
+       -2.00E+01 -3.87E+00 -1.95E+01 -6.16E+00  3.27E+01  3.46E+00 -1.67E-02  0.00E+00  6.74E+00  1.61E+01  3.01E+01
 
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
 #CPUT: Total CPU Time in Seconds,       21.515
Stop Time:
Sat Sep 18 10:13:22 CDT 2021
