Sat Oct 23 23:34:01 CDT 2021
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
$DATA ../../../../data/SD3/SL1/dat66.csv ignore=@
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
Current Date:       23 OCT 2021
Days until program expires : 176
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
 RAW OUTPUT FILE (FILE): m66.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2096.31111326756        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.9378E+02 -4.2372E+01 -1.0181E+01 -1.9241E+01  1.5987E+01  4.3924E+01 -3.2647E+00  1.1084E+01  2.9016E+01 -3.1289E+00
             7.4981E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2104.62901173259        NO. OF FUNC. EVALS.: 159
 CUMULATIVE NO. OF FUNC. EVALS.:      172
 NPARAMETR:  9.4271E-01  1.1598E+00  1.0678E+00  9.8431E-01  1.1202E+00  9.0797E-01  1.0997E+00  8.7347E-01  7.2327E-01  1.0723E+00
             1.0059E+00
 PARAMETER:  4.1000E-02  2.4823E-01  1.6563E-01  8.4184E-02  2.1353E-01  3.4529E-03  1.9505E-01 -3.5282E-02 -2.2397E-01  1.6978E-01
             1.0584E-01
 GRADIENT:  -2.0439E+01  2.6786E+01  4.5245E+00  3.1862E+01  2.1022E+01 -2.4811E+01 -1.6708E+01 -1.3379E+00 -1.6603E+01 -9.6194E+00
             1.3423E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2107.28640826443        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      350
 NPARAMETR:  9.5108E-01  9.2389E-01  1.0650E+00  1.1269E+00  1.0003E+00  9.4317E-01  1.4563E+00  6.6340E-01  5.7776E-01  1.0849E+00
             9.9995E-01
 PARAMETER:  4.9842E-02  2.0841E-02  1.6302E-01  2.1947E-01  1.0026E-01  4.1493E-02  4.7587E-01 -3.1038E-01 -4.4860E-01  1.8151E-01
             9.9949E-02
 GRADIENT:   7.7393E+00  3.1206E+01  5.9948E+00  4.7541E+01  2.6052E+00 -7.6013E+00 -1.0731E+01 -2.1719E+00 -2.1315E+01  1.4044E+00
            -2.1961E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2112.45390290189        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      527
 NPARAMETR:  9.4789E-01  8.2214E-01  8.2744E-01  1.1553E+00  8.1980E-01  9.6021E-01  1.5277E+00  2.6389E-01  7.1831E-01  8.7196E-01
             9.9531E-01
 PARAMETER:  4.6480E-02 -9.5849E-02 -8.9419E-02  2.4440E-01 -9.8691E-02  5.9396E-02  5.2377E-01 -1.2322E+00 -2.3085E-01 -3.7017E-02
             9.5299E-02
 GRADIENT:  -2.3006E+00  8.1453E+00 -7.2895E+00  2.2383E+01  6.6422E+00 -4.4249E-01 -1.9744E+00  3.0278E-01  6.0426E-01  3.7337E-01
            -8.5994E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2112.97511282387        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      702
 NPARAMETR:  9.4696E-01  6.7773E-01  8.6127E-01  1.2268E+00  7.8754E-01  9.5792E-01  1.8036E+00  1.7622E-01  6.8201E-01  8.7152E-01
             1.0002E+00
 PARAMETER:  4.5503E-02 -2.8900E-01 -4.9347E-02  3.0443E-01 -1.3885E-01  5.7004E-02  6.8981E-01 -1.6360E+00 -2.8270E-01 -3.7522E-02
             1.0024E-01
 GRADIENT:   7.4798E-02  4.6825E-01  1.3906E+00 -1.9379E+00 -2.8839E-01 -4.8039E-01  1.3920E-01 -1.0883E-01 -3.6934E-02 -4.2563E-01
            -9.8440E-02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2113.01442348922        NO. OF FUNC. EVALS.: 120
 CUMULATIVE NO. OF FUNC. EVALS.:      822
 NPARAMETR:  9.4716E-01  6.7048E-01  8.6072E-01  1.2293E+00  7.8432E-01  9.5940E-01  1.8182E+00  2.5815E-01  6.7998E-01  8.6271E-01
             9.9966E-01
 PARAMETER:  4.5713E-02 -2.9976E-01 -4.9986E-02  3.0648E-01 -1.4294E-01  5.8553E-02  6.9785E-01 -1.2542E+00 -2.8569E-01 -4.7675E-02
             9.9659E-02
 GRADIENT:   3.6814E+02  4.4258E+01  2.4852E+00  3.1416E+02  1.4046E+01  3.8077E+01  3.6624E+01  2.1023E-01  1.2599E+01  1.0095E+00
             1.2823E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2113.02060257286        NO. OF FUNC. EVALS.: 157
 CUMULATIVE NO. OF FUNC. EVALS.:      979
 NPARAMETR:  9.4620E-01  6.6987E-01  8.6370E-01  1.2313E+00  7.8459E-01  9.5838E-01  1.8171E+00  2.9517E-01  6.7822E-01  8.5826E-01
             9.9913E-01
 PARAMETER:  4.4700E-02 -3.0068E-01 -4.6526E-02  3.0811E-01 -1.4259E-01  5.7489E-02  6.9723E-01 -1.1202E+00 -2.8828E-01 -5.2849E-02
             9.9127E-02
 GRADIENT:  -1.6679E+00  2.1817E-01 -7.4536E-01 -2.0729E+00 -7.1901E-01 -2.7675E-01 -3.8225E-01  1.7835E-02  7.9686E-03  3.2555E-02
             5.3831E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2113.02316595968        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1154
 NPARAMETR:  9.4419E-01  6.6123E-01  8.7234E-01  1.2365E+00  7.8688E-01  9.5712E-01  1.8314E+00  3.1444E-01  6.7507E-01  8.6309E-01
             9.9918E-01
 PARAMETER:  4.2568E-02 -3.1365E-01 -3.6573E-02  3.1232E-01 -1.3968E-01  5.6170E-02  7.0509E-01 -1.0570E+00 -2.9294E-01 -4.7237E-02
             9.9177E-02
 GRADIENT:  -6.3919E+00 -4.0598E-01 -9.4168E-01 -2.9830E+00 -4.4002E-01 -7.5395E-01 -7.9659E-01  3.0575E-02 -1.1813E-01  2.7838E-01
             1.9333E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2113.02514249780        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1329
 NPARAMETR:  9.4290E-01  6.5222E-01  8.8298E-01  1.2417E+00  7.9046E-01  9.5633E-01  1.8504E+00  3.3066E-01  6.7173E-01  8.7015E-01
             9.9935E-01
 PARAMETER:  4.1201E-02 -3.2738E-01 -2.4449E-02  3.1648E-01 -1.3514E-01  5.5351E-02  7.1541E-01 -1.0067E+00 -2.9789E-01 -3.9088E-02
             9.9346E-02
 GRADIENT:  -9.2176E+00 -1.1503E+00 -9.2912E-01 -5.1174E+00  6.0752E-03 -1.0263E+00 -8.6486E-01  2.0577E-02 -1.6577E-01  5.4828E-01
             3.9880E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2113.03198020862        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1507
 NPARAMETR:  9.4085E-01  6.2338E-01  9.2479E-01  1.2609E+00  8.0528E-01  9.5506E-01  1.9175E+00  4.0505E-01  6.6045E-01  8.9202E-01
             9.9939E-01
 PARAMETER:  3.9029E-02 -3.7259E-01  2.1811E-02  3.3180E-01 -1.1657E-01  5.4021E-02  7.5101E-01 -8.0374E-01 -3.1484E-01 -1.4266E-02
             9.9385E-02
 GRADIENT:  -1.2828E+01 -2.0550E+00 -8.6408E-01 -8.0351E+00  1.3655E+00 -1.3635E+00 -9.6964E-01  1.6553E-02 -2.0224E-01  1.0388E+00
             7.3726E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -2113.04264273535        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1682
 NPARAMETR:  9.4078E-01  5.8700E-01  9.8256E-01  1.2876E+00  8.2480E-01  9.5478E-01  2.0123E+00  5.0323E-01  6.4652E-01  9.1499E-01
             9.9897E-01
 PARAMETER:  3.8954E-02 -4.3273E-01  8.2404E-02  3.5281E-01 -9.2609E-02  5.3727E-02  7.9930E-01 -5.8670E-01 -3.3615E-01  1.1156E-02
             9.8967E-02
 GRADIENT:  -1.0917E+01 -1.4670E+00 -5.7275E-01 -6.6249E+00  2.4854E+00 -1.2039E+00 -9.2921E-01  2.8520E-02 -1.5688E-01  1.0718E+00
             7.0055E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -2113.04991753828        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1857
 NPARAMETR:  9.4268E-01  5.5481E-01  1.0292E+00  1.3116E+00  8.3796E-01  9.5581E-01  2.1070E+00  5.7759E-01  6.3591E-01  9.2664E-01
             9.9843E-01
 PARAMETER:  4.0970E-02 -4.8913E-01  1.2875E-01  3.7126E-01 -7.6786E-02  5.4808E-02  8.4527E-01 -4.4889E-01 -3.5271E-01  2.3811E-02
             9.8434E-02
 GRADIENT:  -4.3138E+00 -1.9595E-01  6.1428E-02 -3.0018E+00  2.6806E+00 -5.3057E-01 -5.7998E-01  1.1827E-02 -3.7795E-02  5.6001E-01
             3.2486E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -2113.06481756845        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     2034
 NPARAMETR:  9.4597E-01  5.2433E-01  1.0457E+00  1.3322E+00  8.3449E-01  9.5793E-01  2.2071E+00  6.1157E-01  6.2978E-01  9.1768E-01
             9.9798E-01
 PARAMETER:  4.4455E-02 -5.4564E-01  1.4469E-01  3.8685E-01 -8.0937E-02  5.7016E-02  8.9170E-01 -3.9173E-01 -3.6239E-01  1.4093E-02
             9.7974E-02
 GRADIENT:   5.3652E+00  1.4112E+00  9.4547E-01  2.3713E+00  1.1229E+00  5.1791E-01  8.4025E-02 -4.4538E-02  8.9164E-02 -4.2156E-01
            -3.3882E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -2113.10299993845        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     2211
 NPARAMETR:  9.4451E-01  5.1115E-01  1.0104E+00  1.3357E+00  8.1152E-01  9.5686E-01  2.2444E+00  5.7213E-01  6.3123E-01  8.9750E-01
             9.9877E-01
 PARAMETER:  4.2912E-02 -5.7109E-01  1.1031E-01  3.8942E-01 -1.0885E-01  5.5900E-02  9.0845E-01 -4.5839E-01 -3.6008E-01 -8.1439E-03
             9.8767E-02
 GRADIENT:   1.7595E+00  9.7234E-01 -1.9200E-01  2.9942E+00 -6.6840E-01  1.6948E-01  1.6798E-02  7.0030E-02  1.0132E-02 -1.6071E-03
            -2.7968E-02

0ITERATION NO.:   70    OBJECTIVE VALUE:  -2113.10885235351        NO. OF FUNC. EVALS.: 169
 CUMULATIVE NO. OF FUNC. EVALS.:     2380
 NPARAMETR:  9.4429E-01  5.1216E-01  1.0095E+00  1.3341E+00  8.1202E-01  9.5663E-01  2.2553E+00  5.6862E-01  6.3025E-01  8.9824E-01
             9.9877E-01
 PARAMETER:  4.2677E-02 -5.6911E-01  1.0949E-01  3.8822E-01 -1.0823E-01  5.5664E-02  9.1327E-01 -4.6454E-01 -3.6164E-01 -7.3218E-03
             9.8770E-02
 GRADIENT:  -9.3316E-02  5.0693E-01 -3.2908E-01  4.3884E+00  3.4126E-02 -2.0991E-02  3.8279E-01  2.4007E-02 -4.6372E-02  9.9159E-03
            -3.1047E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     2380
 NO. OF SIG. DIGITS IN FINAL EST.:  2.0

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.2866E-03  2.7973E-02 -2.3836E-02 -3.1132E-02 -6.0943E-03
 SE:             2.9880E-02  2.1063E-02  1.1556E-02  2.2659E-02  2.2501E-02
 N:                     100         100         100         100         100

 P VAL.:         9.6566E-01  1.8416E-01  3.9154E-02  1.6946E-01  7.8651E-01

 ETASHRINKSD(%)  1.0000E-10  2.9437E+01  6.1285E+01  2.4090E+01  2.4619E+01
 ETASHRINKVR(%)  1.0000E-10  5.0208E+01  8.5012E+01  4.2377E+01  4.3178E+01
 EBVSHRINKSD(%)  3.8562E-01  3.0901E+01  6.4204E+01  2.2035E+01  2.2261E+01
 EBVSHRINKVR(%)  7.6976E-01  5.2253E+01  8.7187E+01  3.9214E+01  3.9567E+01
 RELATIVEINF(%)  9.8668E+01  7.8243E+00  1.8102E+00  9.8131E+00  9.6913E+00
 EPSSHRINKSD(%)  3.3309E+01
 EPSSHRINKVR(%)  5.5523E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2113.1088523535145     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1194.1703191488418     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    25.39
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2113.109       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.44E-01  5.12E-01  1.01E+00  1.33E+00  8.12E-01  9.57E-01  2.26E+00  5.69E-01  6.30E-01  8.98E-01  9.99E-01
 


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
 #CPUT: Total CPU Time in Seconds,      192.812
Stop Time:
Sat Oct 23 23:34:29 CDT 2021
