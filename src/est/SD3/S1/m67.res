Sat Oct 23 23:08:48 CDT 2021
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
$DATA ../../../../data/SD3/S1/dat67.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m67.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2109.53593490695        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.3100E+02 -7.8786E+01 -2.7158E+01 -7.0132E+01  3.6319E+01  4.0441E+01 -8.2855E+00  1.3306E+01 -1.2533E+01  1.6921E+01
            -1.1906E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2114.54473957640        NO. OF FUNC. EVALS.: 169
 CUMULATIVE NO. OF FUNC. EVALS.:      182
 NPARAMETR:  9.6568E-01  1.0386E+00  1.0570E+00  1.0024E+00  1.0376E+00  1.0344E+00  1.0068E+00  9.5117E-01  1.0040E+00  9.3116E-01
             1.0091E+00
 PARAMETER:  6.5078E-02  1.3787E-01  1.5548E-01  1.0235E-01  1.3691E-01  1.3385E-01  1.0678E-01  4.9943E-02  1.0400E-01  2.8672E-02
             1.0909E-01
 GRADIENT:  -5.0020E+01 -8.7378E+01 -1.1752E+01 -1.1191E+02  3.3912E+01  1.3232E-01 -1.0406E+01  6.2593E+00 -1.6998E+01 -1.0817E+00
             1.9952E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2118.40169999880        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      357
 NPARAMETR:  9.7795E-01  1.3988E+00  7.4681E-01  8.0750E-01  1.0123E+00  9.4943E-01  1.0862E+00  4.6645E-01  1.2081E+00  8.5567E-01
             1.0315E+00
 PARAMETER:  7.7701E-02  4.3562E-01 -1.9195E-01 -1.1382E-01  1.1223E-01  4.8112E-02  1.8268E-01 -6.6261E-01  2.8906E-01 -5.5870E-02
             1.3105E-01
 GRADIENT:  -3.7084E+01 -3.0086E+00  2.2611E+01 -4.9678E+01 -3.3255E+01 -3.7879E+01  2.2902E+01  1.4617E+00  1.1229E+00  1.7632E+00
             2.2607E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2124.82032591503        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      534
 NPARAMETR:  9.9551E-01  1.3380E+00  6.3259E-01  8.7238E-01  9.0559E-01  1.0380E+00  9.7425E-01  2.7067E-01  1.1713E+00  7.6491E-01
             1.0034E+00
 PARAMETER:  9.5498E-02  3.9120E-01 -3.5794E-01 -3.6533E-02  8.3269E-04  1.3728E-01  7.3916E-02 -1.2069E+00  2.5811E-01 -1.6799E-01
             1.0340E-01
 GRADIENT:   4.2543E+00  9.8664E+00  3.5715E+00  4.0264E+00 -1.5917E+01  3.7911E-01  2.8362E+00  9.9166E-01  4.1612E+00  4.1574E+00
             4.3409E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2125.39608992943        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      710
 NPARAMETR:  9.9361E-01  1.4392E+00  5.7981E-01  8.0091E-01  9.3829E-01  1.0389E+00  9.1149E-01  1.1003E-01  1.2190E+00  7.5506E-01
             9.9645E-01
 PARAMETER:  9.3585E-02  4.6411E-01 -4.4505E-01 -1.2200E-01  3.6302E-02  1.3818E-01  7.3259E-03 -2.1070E+00  2.9807E-01 -1.8096E-01
             9.6443E-02
 GRADIENT:  -5.7743E-01 -1.4021E+00 -1.0105E+00 -8.1238E-01  9.2923E-01  3.9836E-01 -2.3072E-01  1.5768E-01 -1.5700E-01  3.3734E-01
            -1.9753E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2125.47205278573        NO. OF FUNC. EVALS.: 142
 CUMULATIVE NO. OF FUNC. EVALS.:      852
 NPARAMETR:  9.9341E-01  1.4373E+00  5.8091E-01  8.0305E-01  9.3778E-01  1.0375E+00  9.1461E-01  3.0400E-02  1.2197E+00  7.5307E-01
             9.9686E-01
 PARAMETER:  9.3391E-02  4.6277E-01 -4.4317E-01 -1.1934E-01  3.5764E-02  1.3684E-01  1.0744E-02 -3.3933E+00  2.9863E-01 -1.8360E-01
             9.6852E-02
 GRADIENT:  -9.6820E-01 -1.6497E-01 -8.7309E-01  3.8972E-01  1.4274E+00 -1.3763E-01  6.5489E-02  1.2120E-02  2.2045E-01  1.5260E-02
             1.6052E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2125.50831908034        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1028
 NPARAMETR:  9.9438E-01  1.4128E+00  5.8934E-01  8.1757E-01  9.2899E-01  1.0376E+00  9.2310E-01  1.0000E-02  1.2030E+00  7.5018E-01
             9.9737E-01
 PARAMETER:  9.4361E-02  4.4555E-01 -4.2875E-01 -1.0142E-01  2.6337E-02  1.3694E-01  1.9987E-02 -5.9777E+00  2.8482E-01 -1.8744E-01
             9.7362E-02
 GRADIENT:   1.1860E+00 -1.2471E+00 -4.3211E-01 -4.7369E-01  5.4702E-01 -3.4492E-02 -4.7327E-01  0.0000E+00 -9.6869E-02  2.0689E-01
             1.4979E-01

0ITERATION NO.:   34    OBJECTIVE VALUE:  -2125.51524966450        NO. OF FUNC. EVALS.: 127
 CUMULATIVE NO. OF FUNC. EVALS.:     1155
 NPARAMETR:  9.9378E-01  1.4054E+00  5.9235E-01  8.2279E-01  9.2610E-01  1.0377E+00  9.2991E-01  1.0000E-02  1.1981E+00  7.4689E-01
             9.9739E-01
 PARAMETER:  9.3763E-02  4.4030E-01 -4.2367E-01 -9.5058E-02  2.3227E-02  1.3699E-01  2.7333E-02 -6.8033E+00  2.8077E-01 -1.9184E-01
             9.7388E-02
 GRADIENT:   1.7446E-02  2.3617E-02 -1.8730E-03  1.0760E-02 -3.8910E-03 -5.6997E-04 -1.7081E-02  0.0000E+00 -4.7826E-03 -3.7953E-03
            -2.2593E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1155
 NO. OF SIG. DIGITS IN FINAL EST.:  2.2
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY
 THIS MUST BE ADDRESSED BEFORE THE COVARIANCE STEP CAN BE IMPLEMENTED

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.0552E-03 -1.6714E-02 -3.2393E-04  1.1815E-02 -2.2048E-02
 SE:             2.9907E-02  2.3561E-02  1.4067E-04  2.5450E-02  2.0847E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7185E-01  4.7809E-01  2.1289E-02  6.4248E-01  2.9024E-01

 ETASHRINKSD(%)  1.0000E-10  2.1067E+01  9.9529E+01  1.4739E+01  3.0159E+01
 ETASHRINKVR(%)  1.0000E-10  3.7696E+01  9.9998E+01  2.7305E+01  5.1223E+01
 EBVSHRINKSD(%)  3.1869E-01  2.1067E+01  9.9518E+01  1.4693E+01  3.0494E+01
 EBVSHRINKVR(%)  6.3636E-01  3.7695E+01  9.9998E+01  2.7228E+01  5.1689E+01
 RELATIVEINF(%)  9.9295E+01  5.0949E+00  3.9961E-04  7.6214E+00  6.4827E+00
 EPSSHRINKSD(%)  3.3742E+01
 EPSSHRINKVR(%)  5.6099E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2125.5152496644978     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1206.5767164598251     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    12.05
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2125.515       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.94E-01  1.41E+00  5.92E-01  8.23E-01  9.26E-01  1.04E+00  9.30E-01  1.00E-02  1.20E+00  7.47E-01  9.97E-01
 


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
 #CPUT: Total CPU Time in Seconds,       92.123
Stop Time:
Sat Oct 23 23:09:03 CDT 2021
