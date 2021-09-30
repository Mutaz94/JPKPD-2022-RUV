Wed Sep 29 14:12:31 CDT 2021
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
$DATA ../../../../data/spa/S1/dat35.csv ignore=@
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1693.19242683615        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.4468E+02 -2.5017E+00 -5.2993E+01  7.7328E+01  9.8855E+01  6.8786E+01  1.5877E+00  1.1184E+01  2.4897E+00 -1.2651E+01
            -7.4958E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1701.56032513763        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      194
 NPARAMETR:  1.0494E+00  1.0244E+00  1.1296E+00  1.0270E+00  9.8949E-01  9.2976E-01  9.8368E-01  8.8545E-01  1.0382E+00  1.0601E+00
             1.0342E+00
 PARAMETER:  1.4824E-01  1.2412E-01  2.2184E-01  1.2666E-01  8.9438E-02  2.7168E-02  8.3547E-02 -2.1654E-02  1.3745E-01  1.5834E-01
             1.3362E-01
 GRADIENT:  -1.3286E+01  2.9278E+01  1.9796E+00  3.5982E+01 -9.9023E+00 -1.0495E+01  2.3736E+00  3.8716E+00  1.6907E+00 -6.0760E+00
             4.3202E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1702.94694671912        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      372
 NPARAMETR:  1.0499E+00  9.5694E-01  1.2022E+00  1.0530E+00  1.0006E+00  9.4750E-01  8.7342E-01  6.9276E-01  1.0491E+00  1.1631E+00
             1.0205E+00
 PARAMETER:  1.4874E-01  5.5981E-02  2.8412E-01  1.5164E-01  1.0061E-01  4.6072E-02 -3.5339E-02 -2.6707E-01  1.4791E-01  2.5110E-01
             1.2032E-01
 GRADIENT:  -8.7453E+00  1.3057E+01  9.4143E+00  1.2850E+01 -1.2299E+01 -2.3083E+00  1.1161E-01 -1.7484E+00  3.5155E+00 -4.4443E-01
            -2.5275E-02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1703.33295821962        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      551
 NPARAMETR:  1.0514E+00  7.9647E-01  1.2560E+00  1.1530E+00  9.6559E-01  9.4691E-01  1.0088E+00  7.6021E-01  9.5558E-01  1.1337E+00
             1.0206E+00
 PARAMETER:  1.5011E-01 -1.2757E-01  3.2797E-01  2.4236E-01  6.4988E-02  4.5449E-02  1.0875E-01 -1.7416E-01  5.4563E-02  2.2549E-01
             1.2044E-01
 GRADIENT:  -2.1582E+00  6.8674E+00  3.7863E+00  8.8272E+00 -4.0260E+00 -2.0284E+00  2.4997E-01 -3.6047E-01  5.5691E-01 -1.4327E+00
             5.4076E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1703.49865658144        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      727
 NPARAMETR:  1.0509E+00  6.0660E-01  1.3323E+00  1.2742E+00  9.3174E-01  9.5476E-01  1.0500E+00  8.1380E-01  8.8728E-01  1.1508E+00
             1.0151E+00
 PARAMETER:  1.4962E-01 -3.9988E-01  3.8692E-01  3.4230E-01  2.9295E-02  5.3705E-02  1.4879E-01 -1.0604E-01 -1.9597E-02  2.4047E-01
             1.1496E-01
 GRADIENT:   1.8539E+00  5.5631E+00  5.4794E-02  1.1503E+01 -3.2491E+00  2.0583E+00 -5.5192E-01  2.5877E-01 -1.6382E+00  1.2224E+00
            -1.1878E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1703.55908200323        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      902
 NPARAMETR:  1.0490E+00  4.6470E-01  1.4013E+00  1.3640E+00  9.1017E-01  9.5119E-01  1.0920E+00  8.7045E-01  8.4552E-01  1.1575E+00
             1.0157E+00
 PARAMETER:  1.4787E-01 -6.6636E-01  4.3740E-01  4.1045E-01  5.8706E-03  4.9953E-02  1.8799E-01 -3.8742E-02 -6.7799E-02  2.4625E-01
             1.1561E-01
 GRADIENT:   2.1667E+00  5.0329E+00  1.5351E+00  1.1577E+01 -6.8519E+00  1.3310E+00 -4.5444E-01  2.2301E-01 -1.4989E+00  2.1298E+00
            -9.3249E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1703.61069224735        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1079
 NPARAMETR:  1.0471E+00  3.5811E-01  1.4574E+00  1.4307E+00  8.9791E-01  9.4778E-01  1.1393E+00  9.2656E-01  8.1608E-01  1.1563E+00
             1.0170E+00
 PARAMETER:  1.4606E-01 -9.2691E-01  4.7666E-01  4.5815E-01 -7.6900E-03  4.6372E-02  2.3038E-01  2.3728E-02 -1.0325E-01  2.4525E-01
             1.1689E-01
 GRADIENT:   1.4856E+00  3.8065E+00  1.8946E+00  8.9557E+00 -6.9528E+00  5.0039E-01 -2.8807E-01  1.7625E-01 -9.8674E-01  1.9678E+00
            -4.9686E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1703.71927017336        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     1264
 NPARAMETR:  1.0462E+00  3.2296E-01  1.4694E+00  1.4471E+00  8.9696E-01  9.4641E-01  1.3671E+00  9.3747E-01  8.0547E-01  1.1450E+00
             1.0177E+00
 PARAMETER:  1.4517E-01 -1.0302E+00  4.8488E-01  4.6959E-01 -8.7493E-03  4.4924E-02  4.1267E-01  3.5432E-02 -1.1632E-01  2.3542E-01
             1.1752E-01
 GRADIENT:   6.9436E-01  8.5437E-01 -1.6727E-01 -4.7898E+00  8.0447E-01  8.9429E-02 -4.5186E-02 -6.1188E-02  1.2200E+00  3.1314E-01
            -1.0430E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1703.73767611532        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1447             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0469E+00  3.0729E-01  1.4718E+00  1.4535E+00  8.9374E-01  9.4619E-01  1.5608E+00  9.4453E-01  7.9674E-01  1.1394E+00
             1.0182E+00
 PARAMETER:  1.4581E-01 -1.0800E+00  4.8649E-01  4.7397E-01 -1.2338E-02  4.4684E-02  5.4519E-01  4.2929E-02 -1.2722E-01  2.3054E-01
             1.1800E-01
 GRADIENT:   6.4770E+02  4.4581E+01  7.6204E+00  7.4708E+02  9.3428E+00  4.6143E+01  3.7552E+00  2.1028E-01  1.1541E+01  2.5278E+00
             1.1529E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1703.74326074716        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1631             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0472E+00  3.0660E-01  1.4702E+00  1.4545E+00  8.9227E-01  9.4634E-01  1.5774E+00  9.4434E-01  7.9448E-01  1.1376E+00
             1.0181E+00
 PARAMETER:  1.4607E-01 -1.0822E+00  4.8542E-01  4.7463E-01 -1.3984E-02  4.4842E-02  5.5581E-01  4.2732E-02 -1.3006E-01  2.2891E-01
             1.1793E-01
 GRADIENT:   6.4952E+02  4.4793E+01  7.9702E+00  7.5018E+02  8.5000E+00  4.6169E+01  3.8422E+00  2.2401E-01  1.1107E+01  2.4840E+00
             1.0457E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1703.74620578813        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1808             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0471E+00  3.0540E-01  1.4671E+00  1.4553E+00  8.9100E-01  9.4633E-01  1.5769E+00  9.4318E-01  7.9332E-01  1.1358E+00
             1.0181E+00
 PARAMETER:  1.4605E-01 -1.0861E+00  4.8332E-01  4.7524E-01 -1.5416E-02  4.4839E-02  5.5548E-01  4.1497E-02 -1.3153E-01  2.2736E-01
             1.1795E-01
 GRADIENT:   6.4919E+02  4.4691E+01  7.7316E+00  7.5243E+02  8.6528E+00  4.6155E+01  3.7681E+00  2.6530E-01  1.0831E+01  2.3754E+00
             1.0132E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1703.74817805075        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1990             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0471E+00  3.0412E-01  1.4641E+00  1.4561E+00  8.8977E-01  9.4633E-01  1.5804E+00  9.4177E-01  7.9228E-01  1.1344E+00
             1.0181E+00
 PARAMETER:  1.4603E-01 -1.0903E+00  4.8126E-01  4.7575E-01 -1.6798E-02  4.4838E-02  5.5768E-01  4.0001E-02 -1.3283E-01  2.2612E-01
             1.1792E-01
 GRADIENT:   6.4897E+02  4.4537E+01  7.5017E+00  7.5424E+02  8.8397E+00  4.6145E+01  3.7388E+00  3.0363E-01  1.0642E+01  2.3179E+00
             9.8422E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1703.74974246709        NO. OF FUNC. EVALS.: 174
 CUMULATIVE NO. OF FUNC. EVALS.:     2164
 NPARAMETR:  1.0466E+00  3.0277E-01  1.4640E+00  1.4571E+00  8.8849E-01  9.4619E-01  1.5910E+00  9.3808E-01  7.9217E-01  1.1339E+00
             1.0180E+00
 PARAMETER:  1.4553E-01 -1.0948E+00  4.8119E-01  4.7645E-01 -1.8230E-02  4.4688E-02  5.6433E-01  3.6081E-02 -1.3298E-01  2.2566E-01
             1.1788E-01
 GRADIENT:   2.3421E+00  4.2663E-01  4.8191E-01 -1.0232E+01  2.4045E-01  6.9742E-02  7.0554E-03 -2.3969E-02  1.4388E-01  7.0024E-02
            -2.2813E-02

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1703.75175874990        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     2341
 NPARAMETR:  1.0471E+00  3.0113E-01  1.4623E+00  1.4575E+00  8.8716E-01  9.4632E-01  1.6160E+00  9.3679E-01  7.9160E-01  1.1331E+00
             1.0180E+00
 PARAMETER:  1.4600E-01 -1.1002E+00  4.7998E-01  4.7672E-01 -1.9731E-02  4.4824E-02  5.7996E-01  3.4703E-02 -1.3370E-01  2.2498E-01
             1.1788E-01
 GRADIENT:   3.5798E+00  3.0971E-01  7.0957E-01 -1.1410E+01 -1.5165E-01  1.2547E-01  3.3484E-02 -1.5483E-02  3.1306E-01  1.7302E-01
            -1.9244E-03

0ITERATION NO.:   67    OBJECTIVE VALUE:  -1703.75199072224        NO. OF FUNC. EVALS.:  59
 CUMULATIVE NO. OF FUNC. EVALS.:     2400
 NPARAMETR:  1.0467E+00  3.0079E-01  1.4619E+00  1.4581E+00  8.8719E-01  9.4622E-01  1.6067E+00  9.3689E-01  7.9147E-01  1.1329E+00
             1.0180E+00
 PARAMETER:  1.4564E-01 -1.1014E+00  4.7974E-01  4.7713E-01 -1.9697E-02  4.4715E-02  5.7416E-01  3.4815E-02 -1.3387E-01  2.2482E-01
             1.1788E-01
 GRADIENT:  -8.4678E-01  1.6548E-01  4.8032E-01  3.9206E-01  1.6697E-01 -4.1861E-02 -5.1919E-03 -3.0706E-03  1.5551E-01  1.0574E-01
            -1.3083E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     2400
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.3253E-04 -7.6616E-03 -2.3841E-02 -4.2912E-03 -3.0940E-02
 SE:             2.9854E-02  7.8819E-03  1.3176E-02  2.8291E-02  2.3290E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9379E-01  3.3103E-01  7.0379E-02  8.7944E-01  1.8402E-01

 ETASHRINKSD(%)  1.0000E-10  7.3594E+01  5.5859E+01  5.2225E+00  2.1976E+01
 ETASHRINKVR(%)  1.0000E-10  9.3027E+01  8.0516E+01  1.0172E+01  3.9122E+01
 EBVSHRINKSD(%)  4.6550E-01  7.3922E+01  5.9497E+01  5.3618E+00  1.8204E+01
 EBVSHRINKVR(%)  9.2884E-01  9.3200E+01  8.3595E+01  1.0436E+01  3.3095E+01
 RELATIVEINF(%)  9.5529E+01  1.6812E-01  3.5703E+00  2.6437E+00  7.0435E+00
 EPSSHRINKSD(%)  4.4027E+01
 EPSSHRINKVR(%)  6.8670E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1703.7519907222427     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -968.60116415850450     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    31.73
 Elapsed covariance  time in seconds:     6.06
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1703.752       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.05E+00  3.01E-01  1.46E+00  1.46E+00  8.87E-01  9.46E-01  1.61E+00  9.37E-01  7.91E-01  1.13E+00  1.02E+00
 


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
 
         3.55E-02  8.35E-01  5.39E-01  5.04E-01  2.13E-01  6.27E-02  4.18E+00  5.72E-01  3.18E-01  2.37E-01  7.65E-02
 


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
+        1.26E-03
 
 TH 2
+        1.67E-02  6.97E-01
 
 TH 3
+       -7.13E-03 -2.42E-01  2.91E-01
 
 TH 4
+       -1.01E-02 -4.20E-01  1.58E-01  2.54E-01
 
 TH 5
+        2.20E-03  1.13E-01  3.29E-02 -6.37E-02  4.52E-02
 
 TH 6
+        3.22E-04  6.22E-03 -2.42E-03 -3.99E-03  1.28E-03  3.93E-03
 
 TH 7
+       -7.27E-02 -3.00E+00  4.96E-01  1.77E+00 -6.97E-01 -2.92E-02  1.75E+01
 
 TH 8
+       -8.66E-03 -3.24E-01  2.73E-01  2.06E-01  1.31E-03 -4.50E-03  1.02E+00  3.27E-01
 
 TH 9
+        5.58E-03  2.49E-01 -7.17E-02 -1.49E-01  4.66E-02  2.32E-03 -1.25E+00 -1.06E-01  1.01E-01
 
 TH10
+        8.60E-04  6.20E-02  6.40E-02 -3.20E-02  4.23E-02  1.42E-03 -5.43E-01  2.60E-02  3.12E-02  5.63E-02
 
 TH11
+        1.78E-04  6.13E-03  2.60E-03 -3.42E-03  2.36E-03 -6.41E-04 -5.71E-02  5.90E-04  3.66E-03  2.01E-04  5.85E-03
 
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
+        3.55E-02
 
 TH 2
+        5.63E-01  8.35E-01
 
 TH 3
+       -3.72E-01 -5.37E-01  5.39E-01
 
 TH 4
+       -5.62E-01 -9.96E-01  5.81E-01  5.04E-01
 
 TH 5
+        2.91E-01  6.38E-01  2.87E-01 -5.94E-01  2.13E-01
 
 TH 6
+        1.45E-01  1.19E-01 -7.15E-02 -1.26E-01  9.61E-02  6.27E-02
 
 TH 7
+       -4.90E-01 -8.60E-01  2.20E-01  8.42E-01 -7.84E-01 -1.11E-01  4.18E+00
 
 TH 8
+       -4.27E-01 -6.78E-01  8.86E-01  7.13E-01  1.08E-02 -1.25E-01  4.26E-01  5.72E-01
 
 TH 9
+        4.94E-01  9.40E-01 -4.19E-01 -9.32E-01  6.89E-01  1.16E-01 -9.39E-01 -5.86E-01  3.18E-01
 
 TH10
+        1.02E-01  3.13E-01  5.01E-01 -2.68E-01  8.38E-01  9.56E-02 -5.48E-01  1.91E-01  4.13E-01  2.37E-01
 
 TH11
+        6.54E-02  9.59E-02  6.29E-02 -8.86E-02  1.45E-01 -1.34E-01 -1.79E-01  1.35E-02  1.50E-01  1.11E-02  7.65E-02
 
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
+        1.28E+03
 
 TH 2
+       -8.96E+01  4.01E+02
 
 TH 3
+        4.65E+00  6.08E+01  1.20E+02
 
 TH 4
+       -5.96E+01  5.68E+02 -4.56E+00  9.57E+02
 
 TH 5
+        7.07E+01 -2.48E+02 -2.44E+02 -8.54E+01  8.53E+02
 
 TH 6
+       -6.82E+01  5.11E+01  1.17E+01  6.18E+01 -5.40E+01  2.77E+02
 
 TH 7
+        7.43E+00  5.54E-02  1.21E-01  1.54E+00  3.50E+00 -2.00E-01  7.82E-01
 
 TH 8
+       -4.66E+00 -1.33E+01 -2.87E+01 -2.16E+01  7.56E+00  4.09E+00 -1.76E-01  2.45E+01
 
 TH 9
+        1.24E+02 -3.78E+00  1.82E+01  5.44E+01  2.18E+00 -2.58E+00  9.58E+00 -5.40E+00  2.24E+02
 
 TH10
+       -6.39E+00  7.55E+00 -1.89E+01 -1.75E+01 -1.10E+02 -2.26E+00  2.62E-01  1.92E+01 -1.80E+01  1.07E+02
 
 TH11
+       -2.22E+01 -2.93E+00 -2.68E+01 -6.83E+00  1.99E+00  3.09E+01  7.94E-01  1.01E+01 -2.24E+01  4.26E+01  2.04E+02
 
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
 #CPUT: Total CPU Time in Seconds,       37.840
Stop Time:
Wed Sep 29 14:13:10 CDT 2021
