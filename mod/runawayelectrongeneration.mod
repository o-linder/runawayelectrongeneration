	  G  º   k820309    j          14.0        Ö¡l_                                                                                                           
       src/runawayelectrongeneration.f RUNAWAYELECTRONGENERATION              HOT_TAIL_DENSITY DREICER_GROWTHRATE_CLASSIC DREICER_GROWTHRATE_CODE_NEURAL_NETWORK AVALANCHE_GROWTHRATE_CLASSIC AVALANCHE_GROWTHRATE E_CEFF_OVER_E_C E_C E_D TAU EVDF_MODIFIED_MAXWELLIAN P_STAR NUBAR_D CALC_NUBAR_D01 NUBAR_S CALC_NUBAR_S01 NU_EE V_C V_TH LN_LAMBDA_0 LN_LAMBDA_C LN_LAMBDA_EE LN_LAMBDA_EI                      @                              
       HOT_TAIL_DENSITY EVDF_MODIFIED_MAXWELLIAN TAU V_C V_TH                      @                              
       DREICER_GROWTHRATE_CLASSIC DREICER_GROWTHRATE_CODE_NEURAL_NETWORK                      @                              
       AVALANCHE_GROWTHRATE_CLASSIC AVALANCHE_GROWTHRATE P_STAR                      @                              
       E_CEFF_OVER_E_C E_C E_D                      @                              
       NU_EE NUBAR_D CALC_NUBAR_D01 NUBAR_S CALC_NUBAR_S01                      @                              
       LN_LAMBDA_0 LN_LAMBDA_C LN_LAMBDA_EE LN_LAMBDA_EI %         @                                                  
       #HOT_TAIL_DENSITY%REAL    #HOT_TAIL_DENSITY%PRESENT 	   #HOT_TAIL_DENSITY%SQRT 
   #HOT_TAIL_DENSITY%LOG    #HOT_TAIL_DENSITY%ABS    #HOT_TAIL_DENSITY%MOD    #HOT_TAIL_DENSITY%EXP    #TIME    #T_DEC    #EPAR    #T    #T_I    #T_F    #NE    #NE_I    #NE_F    #VERBOSE    #STORE_RESULT    #UN                  @              @                  REAL               @                            	     PRESENT               @                            
     SQRT               @                                 LOG               @                                 ABS               @                                 MOD               @                                 EXP           
                                      
                
                                      
                
                                      
                
                                     
                
                                     
                
                                     
                
                                     
                
                                     
                
                                     
                
                                                      
                                                      
                                            %         @                                                  
       #EVDF_MODIFIED_MAXWELLIAN%SQRT    #EVDF_MODIFIED_MAXWELLIAN%EXP    #V    #NE    #V_TH     #TAU !                 @                                 SQRT               @                                 EXP           
                                      
                
                                      
                
                                       
                
                                 !     
      %         @                                "                    
       #TIME #   #T_CHAR $   #NE_I %   #NE_F &   #T_I '             
                                 #     
                
                                 $     
                
                                 %     
                
                                 &     
                
                                 '     
      %         @                                (                   
       #V_C%SQRT )   #NE *   #T +   #EPAR ,                 @                            )     SQRT           
                                 *     
                
                                 +     
                
                                 ,     
      %         @                                -                   
       #V_TH%SQRT .   #T /                 @                            .     SQRT           
                                 /     
      %         @                               0                   
       #DREICER_GROWTHRATE_CLASSIC%SUM 1   #DREICER_GROWTHRATE_CLASSIC%SQRT 2   #DREICER_GROWTHRATE_CLASSIC%ASIN 3   #DREICER_GROWTHRATE_CLASSIC%EXP 4   #Z 5   #Z0 6   #NJ 7   #NSPECIES 8   #T 9   #EPAR :                 @                            1     SUM               @                            2     SQRT               @                            3     ASIN               @                            4     EXP          
                                  5                        p          5 O p            5 O p                                   
                                  6                        p          5 O p            5 O p                                   
                                 7                    
    p          5 O p            5 O p                                    
                                  8                     
                                 9     
                
                                 :     
      %         @                               ;                   
       #DREICER_GROWTHRATE_CODE_NEURAL_NETWORK%SUM <   #DREICER_GROWTHRATE_CODE_NEURAL_NETWORK%SQRT =   #DREICER_GROWTHRATE_CODE_NEURAL_NETWORK%EXP >   #DREICER_GROWTHRATE_CODE_NEURAL_NETWORK%RESHAPE ?   #DREICER_GROWTHRATE_CODE_NEURAL_NETWORK%LOG @   #DREICER_GROWTHRATE_CODE_NEURAL_NETWORK%DOT_PRODUCT A   #DREICER_GROWTHRATE_CODE_NEURAL_NETWORK%TANH B   #DREICER_GROWTHRATE_CODE_NEURAL_NETWORK%MATMUL C   #Z D   #Z0 E   #NJ F   #NSPECIES G   #T H   #EPAR I                 @                            <     SUM               @                            =     SQRT               @                            >     EXP               @                            ?     RESHAPE               @                            @     LOG               @                            A     DOT_PRODUCT               @                            B     TANH               @                            C     MATMUL          
                                  D                        p          5 O p            5 O p                                   
                                  E                        p          5 O p            5 O p                                   
                                 F                    
    p          5 O p            5 O p                                    
                                  G                     
                                 H     
                
                                 I     
      %         @                               J                   
       #AVALANCHE_GROWTHRATE_CLASSIC%SUM K   #AVALANCHE_GROWTHRATE_CLASSIC%SQRT L   #Z M   #Z0 N   #NJ O   #NSPECIES P   #T Q   #EPAR R   #EPS S                 @                            K     SUM               @                            L     SQRT          
                                  M                        p          5 O p            5 O p                                   
                                  N                        p          5 O p            5 O p                                   
                                 O                    
    p          5 O p            5 O p                                    
                                  P                     
                                 Q     
                
                                 R     
                
                                 S     
      %         @                               T                   
       #AVALANCHE_GROWTHRATE%SUM U   #AVALANCHE_GROWTHRATE%SQRT V   #Z W   #Z0 X   #NJ Y   #NSPECIES Z   #T [   #EPAR \   #B ]                 @                            U     SUM               @                            V     SQRT          
                                  W                        p          5 O p            5 O p                                   
                                  X                        p          5 O p            5 O p                                   
                                 Y                    
    p          5 O p            5 O p                                    
                                  Z                     
                                 [     
                
                                 \     
                
                                 ]     
      %         @                                ^                   
       #P_STAR%SUM _   #P_STAR%SQRT `   #P_STAR%MAX a   #P_STAR%ABS b   #Z c   #Z0 d   #NJ e   #NSPECIES f   #T g   #EPAR h   #B i                 @                            _     SUM               @                            `     SQRT               @                            a     MAX               @                            b     ABS          
                                  c                        p          5 O p            5 O p                                   
                                  d                        p          5 O p            5 O p                                   
                                 e                    
 	   p          5 O p            5 O p                                    
                                  f                     
                                 g     
                
                                 h     
                
                                 i     
      %         @                                j                   
       #E_CEFF_OVER_E_C%SUM k   #E_CEFF_OVER_E_C%LOG l   #E_CEFF_OVER_E_C%SQRT m   #E_CEFF_OVER_E_C%ABS n   #Z o   #Z0 p   #NJ q   #NSPECIES r   #T s   #B t                 @                            k     SUM               @                            l     LOG               @                            m     SQRT               @                            n     ABS          
                                  o                        p          5 O p            5 O p                                   
                                  p                        p          5 O p            5 O p                                   
                                 q                    
    p          5 O p            5 O p                                    
                                  r                     
                                 s     
                
                                 t     
      %         @                                u                    
       #NE v   #T w             
                                 v     
                
                                 w     
      %         @                                x                    
       #NE y   #T z             
                                 y     
                
                                 z     
      %         @                                {                    
       #NE |   #T }             
                                 |     
                
                                 }     
      %         @                                ~                   
       #NUBAR_D%SUM    #NUBAR_D%LOG    #Z    #Z0    #NJ    #NSPECIES    #T    #P                  @                                 SUM               @                                 LOG          
                                                          p          5 O p            5 O p                                   
                                                          p          5 O p            5 O p                                   
                                                     
    p          5 O p            5 O p                                    
                                                       
                                      
                
                                      
      #         @                                                      #CALC_NUBAR_D01%SUM    #CALC_NUBAR_D01%LOG    #Z    #Z0    #NJ    #NSPECIES    #T    #NUBAR_D0    #NUBAR_D1                  @                                 SUM               @                                 LOG          
                                                          p          5 O p            5 O p                                   
                                                          p          5 O p            5 O p                                   
                                                     
    p          5 O p            5 O p                                    
                                                       
                                      
                                                     
                                                      
       %         @                                                   
       #NUBAR_S%SUM    #NUBAR_S%LOG    #NUBAR_S%SQRT    #Z    #Z0    #NJ    #NSPECIES    #T    #P                  @                                 SUM               @                                 LOG               @                                 SQRT          
                                                          p          5 O p            5 O p                                   
                                                          p          5 O p            5 O p                                   
                                                     
 	   p          5 O p            5 O p                                    
                                                       
                                      
                
                                      
      #         @                                                      #CALC_NUBAR_S01%SUM    #CALC_NUBAR_S01%LOG    #Z    #Z0    #NJ     #NSPECIES ¡   #T ¢   #NUBAR_S0 £   #NUBAR_S1 ¤                 @                                 SUM               @                                 LOG          
                                                       
   p          5 O p            5 O p                                   
                                                          p          5 O p            5 O p                                   
                                                      
    p          5 O p            5 O p                                    
                                  ¡                     
                                 ¢     
                                                £     
                                                 ¤     
       %         @                                ¥                   
       #LN_LAMBDA_0%LOG ¦   #NE §   #T ¨                 @                            ¦     LOG           
                                 §     
                
                                 ¨     
      %         @                                ©                   
       #LN_LAMBDA_C%LOG ª   #NE «   #T ¬                 @                            ª     LOG           
                                 «     
                
                                 ¬     
      %         @                                ­                   
       #LN_LAMBDA_EE%LOG ®   #LN_LAMBDA_EE%SQRT ¯   #NE °   #T ±   #P ²                 @                            ®     LOG               @                            ¯     SQRT           
                                 °     
                
                                 ±     
                
                                 ²     
      %         @                                ³                   
       #LN_LAMBDA_EI%LOG ´   #LN_LAMBDA_EI%SQRT µ   #NE ¶   #T ·   #P ¸                 @                            ´     LOG               @                            µ     SQRT           
                                 ¶     
                
                                 ·     
                
                                 ¸     
             B      fn#fn /   â   @  b   uapp(RUNAWAYELECTRONGENERATION .   "  w   J  CALCULATE_HOT_TAIL_POPULATION -        J  CALCULATE_DREICER_GROWTHRATE /     y   J  CALCULATE_AVALANCHE_GROWTHRATE       X   J  ELECTRIC_FIELDS &   ì  t   J  COLLISION_FREQUENCIES #   `  r   J  COULOMB_LOGARITHMS ?   Ò        HOT_TAIL_DENSITY+CALCULATE_HOT_TAIL_POPULATION I   Y  =      HOT_TAIL_DENSITY%REAL+CALCULATE_HOT_TAIL_POPULATION=REAL O     @      HOT_TAIL_DENSITY%PRESENT+CALCULATE_HOT_TAIL_POPULATION=PRESENT I   Ö  =      HOT_TAIL_DENSITY%SQRT+CALCULATE_HOT_TAIL_POPULATION=SQRT G     <      HOT_TAIL_DENSITY%LOG+CALCULATE_HOT_TAIL_POPULATION=LOG G   O  <      HOT_TAIL_DENSITY%ABS+CALCULATE_HOT_TAIL_POPULATION=ABS G     <      HOT_TAIL_DENSITY%MOD+CALCULATE_HOT_TAIL_POPULATION=MOD G   Ç  <      HOT_TAIL_DENSITY%EXP+CALCULATE_HOT_TAIL_POPULATION=EXP D     @   a   HOT_TAIL_DENSITY%TIME+CALCULATE_HOT_TAIL_POPULATION E   C  @   a   HOT_TAIL_DENSITY%T_DEC+CALCULATE_HOT_TAIL_POPULATION D     @   a   HOT_TAIL_DENSITY%EPAR+CALCULATE_HOT_TAIL_POPULATION A   Ã  @   a   HOT_TAIL_DENSITY%T+CALCULATE_HOT_TAIL_POPULATION C   	  @   a   HOT_TAIL_DENSITY%T_I+CALCULATE_HOT_TAIL_POPULATION C   C	  @   a   HOT_TAIL_DENSITY%T_F+CALCULATE_HOT_TAIL_POPULATION B   	  @   a   HOT_TAIL_DENSITY%NE+CALCULATE_HOT_TAIL_POPULATION D   Ã	  @   a   HOT_TAIL_DENSITY%NE_I+CALCULATE_HOT_TAIL_POPULATION D   
  @   a   HOT_TAIL_DENSITY%NE_F+CALCULATE_HOT_TAIL_POPULATION G   C
  @   a   HOT_TAIL_DENSITY%VERBOSE+CALCULATE_HOT_TAIL_POPULATION L   
  @   a   HOT_TAIL_DENSITY%STORE_RESULT+CALCULATE_HOT_TAIL_POPULATION B   Ã
  @   a   HOT_TAIL_DENSITY%UN+CALCULATE_HOT_TAIL_POPULATION G     ·       EVDF_MODIFIED_MAXWELLIAN+CALCULATE_HOT_TAIL_POPULATION Q   º  =      EVDF_MODIFIED_MAXWELLIAN%SQRT+CALCULATE_HOT_TAIL_POPULATION=SQRT O   ÷  <      EVDF_MODIFIED_MAXWELLIAN%EXP+CALCULATE_HOT_TAIL_POPULATION=EXP I   3  @   a   EVDF_MODIFIED_MAXWELLIAN%V+CALCULATE_HOT_TAIL_POPULATION J   s  @   a   EVDF_MODIFIED_MAXWELLIAN%NE+CALCULATE_HOT_TAIL_POPULATION L   ³  @   a   EVDF_MODIFIED_MAXWELLIAN%V_TH+CALCULATE_HOT_TAIL_POPULATION K   ó  @   a   EVDF_MODIFIED_MAXWELLIAN%TAU+CALCULATE_HOT_TAIL_POPULATION 2   3         TAU+CALCULATE_HOT_TAIL_POPULATION 7   ¶  @   a   TAU%TIME+CALCULATE_HOT_TAIL_POPULATION 9   ö  @   a   TAU%T_CHAR+CALCULATE_HOT_TAIL_POPULATION 7   6  @   a   TAU%NE_I+CALCULATE_HOT_TAIL_POPULATION 7   v  @   a   TAU%NE_F+CALCULATE_HOT_TAIL_POPULATION 6   ¶  @   a   TAU%T_I+CALCULATE_HOT_TAIL_POPULATION 2   ö  w       V_C+CALCULATE_HOT_TAIL_POPULATION <   m  =      V_C%SQRT+CALCULATE_HOT_TAIL_POPULATION=SQRT 5   ª  @   a   V_C%NE+CALCULATE_HOT_TAIL_POPULATION 4   ê  @   a   V_C%T+CALCULATE_HOT_TAIL_POPULATION 7   *  @   a   V_C%EPAR+CALCULATE_HOT_TAIL_POPULATION 3   j  f       V_TH+CALCULATE_HOT_TAIL_POPULATION =   Ð  =      V_TH%SQRT+CALCULATE_HOT_TAIL_POPULATION=SQRT 5     @   a   V_TH%T+CALCULATE_HOT_TAIL_POPULATION H   M        DREICER_GROWTHRATE_CLASSIC+CALCULATE_DREICER_GROWTHRATE P   e  <      DREICER_GROWTHRATE_CLASSIC%SUM+CALCULATE_DREICER_GROWTHRATE=SUM R   ¡  =      DREICER_GROWTHRATE_CLASSIC%SQRT+CALCULATE_DREICER_GROWTHRATE=SQRT R   Þ  =      DREICER_GROWTHRATE_CLASSIC%ASIN+CALCULATE_DREICER_GROWTHRATE=ASIN P     <      DREICER_GROWTHRATE_CLASSIC%EXP+CALCULATE_DREICER_GROWTHRATE=EXP J   W  ¤   a   DREICER_GROWTHRATE_CLASSIC%Z+CALCULATE_DREICER_GROWTHRATE K   û  ¤   a   DREICER_GROWTHRATE_CLASSIC%Z0+CALCULATE_DREICER_GROWTHRATE K     ¤   a   DREICER_GROWTHRATE_CLASSIC%NJ+CALCULATE_DREICER_GROWTHRATE Q   C  @   a   DREICER_GROWTHRATE_CLASSIC%NSPECIES+CALCULATE_DREICER_GROWTHRATE J     @   a   DREICER_GROWTHRATE_CLASSIC%T+CALCULATE_DREICER_GROWTHRATE M   Ã  @   a   DREICER_GROWTHRATE_CLASSIC%EPAR+CALCULATE_DREICER_GROWTHRATE T           DREICER_GROWTHRATE_CODE_NEURAL_NETWORK+CALCULATE_DREICER_GROWTHRATE \     <      DREICER_GROWTHRATE_CODE_NEURAL_NETWORK%SUM+CALCULATE_DREICER_GROWTHRATE=SUM ^   V  =      DREICER_GROWTHRATE_CODE_NEURAL_NETWORK%SQRT+CALCULATE_DREICER_GROWTHRATE=SQRT \     <      DREICER_GROWTHRATE_CODE_NEURAL_NETWORK%EXP+CALCULATE_DREICER_GROWTHRATE=EXP d   Ï  @      DREICER_GROWTHRATE_CODE_NEURAL_NETWORK%RESHAPE+CALCULATE_DREICER_GROWTHRATE=RESHAPE \     <      DREICER_GROWTHRATE_CODE_NEURAL_NETWORK%LOG+CALCULATE_DREICER_GROWTHRATE=LOG l   K  D      DREICER_GROWTHRATE_CODE_NEURAL_NETWORK%DOT_PRODUCT+CALCULATE_DREICER_GROWTHRATE=DOT_PRODUCT ^     =      DREICER_GROWTHRATE_CODE_NEURAL_NETWORK%TANH+CALCULATE_DREICER_GROWTHRATE=TANH b   Ì  ?      DREICER_GROWTHRATE_CODE_NEURAL_NETWORK%MATMUL+CALCULATE_DREICER_GROWTHRATE=MATMUL V     ¤   a   DREICER_GROWTHRATE_CODE_NEURAL_NETWORK%Z+CALCULATE_DREICER_GROWTHRATE W   ¯  ¤   a   DREICER_GROWTHRATE_CODE_NEURAL_NETWORK%Z0+CALCULATE_DREICER_GROWTHRATE W   S  ¤   a   DREICER_GROWTHRATE_CODE_NEURAL_NETWORK%NJ+CALCULATE_DREICER_GROWTHRATE ]   ÷  @   a   DREICER_GROWTHRATE_CODE_NEURAL_NETWORK%NSPECIES+CALCULATE_DREICER_GROWTHRATE V   7  @   a   DREICER_GROWTHRATE_CODE_NEURAL_NETWORK%T+CALCULATE_DREICER_GROWTHRATE Y   w  @   a   DREICER_GROWTHRATE_CODE_NEURAL_NETWORK%EPAR+CALCULATE_DREICER_GROWTHRATE L   ·  Ü       AVALANCHE_GROWTHRATE_CLASSIC+CALCULATE_AVALANCHE_GROWTHRATE T     <      AVALANCHE_GROWTHRATE_CLASSIC%SUM+CALCULATE_AVALANCHE_GROWTHRATE=SUM V   Ï  =      AVALANCHE_GROWTHRATE_CLASSIC%SQRT+CALCULATE_AVALANCHE_GROWTHRATE=SQRT N     ¤   a   AVALANCHE_GROWTHRATE_CLASSIC%Z+CALCULATE_AVALANCHE_GROWTHRATE O   °  ¤   a   AVALANCHE_GROWTHRATE_CLASSIC%Z0+CALCULATE_AVALANCHE_GROWTHRATE O   T  ¤   a   AVALANCHE_GROWTHRATE_CLASSIC%NJ+CALCULATE_AVALANCHE_GROWTHRATE U   ø  @   a   AVALANCHE_GROWTHRATE_CLASSIC%NSPECIES+CALCULATE_AVALANCHE_GROWTHRATE N   8   @   a   AVALANCHE_GROWTHRATE_CLASSIC%T+CALCULATE_AVALANCHE_GROWTHRATE Q   x   @   a   AVALANCHE_GROWTHRATE_CLASSIC%EPAR+CALCULATE_AVALANCHE_GROWTHRATE P   ¸   @   a   AVALANCHE_GROWTHRATE_CLASSIC%EPS+CALCULATE_AVALANCHE_GROWTHRATE D   ø   Ê       AVALANCHE_GROWTHRATE+CALCULATE_AVALANCHE_GROWTHRATE L   Â!  <      AVALANCHE_GROWTHRATE%SUM+CALCULATE_AVALANCHE_GROWTHRATE=SUM N   þ!  =      AVALANCHE_GROWTHRATE%SQRT+CALCULATE_AVALANCHE_GROWTHRATE=SQRT F   ;"  ¤   a   AVALANCHE_GROWTHRATE%Z+CALCULATE_AVALANCHE_GROWTHRATE G   ß"  ¤   a   AVALANCHE_GROWTHRATE%Z0+CALCULATE_AVALANCHE_GROWTHRATE G   #  ¤   a   AVALANCHE_GROWTHRATE%NJ+CALCULATE_AVALANCHE_GROWTHRATE M   '$  @   a   AVALANCHE_GROWTHRATE%NSPECIES+CALCULATE_AVALANCHE_GROWTHRATE F   g$  @   a   AVALANCHE_GROWTHRATE%T+CALCULATE_AVALANCHE_GROWTHRATE I   §$  @   a   AVALANCHE_GROWTHRATE%EPAR+CALCULATE_AVALANCHE_GROWTHRATE F   ç$  @   a   AVALANCHE_GROWTHRATE%B+CALCULATE_AVALANCHE_GROWTHRATE 6   '%  Î       P_STAR+CALCULATE_AVALANCHE_GROWTHRATE >   õ%  <      P_STAR%SUM+CALCULATE_AVALANCHE_GROWTHRATE=SUM @   1&  =      P_STAR%SQRT+CALCULATE_AVALANCHE_GROWTHRATE=SQRT >   n&  <      P_STAR%MAX+CALCULATE_AVALANCHE_GROWTHRATE=MAX >   ª&  <      P_STAR%ABS+CALCULATE_AVALANCHE_GROWTHRATE=ABS 8   æ&  ¤   a   P_STAR%Z+CALCULATE_AVALANCHE_GROWTHRATE 9   '  ¤   a   P_STAR%Z0+CALCULATE_AVALANCHE_GROWTHRATE 9   .(  ¤   a   P_STAR%NJ+CALCULATE_AVALANCHE_GROWTHRATE ?   Ò(  @   a   P_STAR%NSPECIES+CALCULATE_AVALANCHE_GROWTHRATE 8   )  @   a   P_STAR%T+CALCULATE_AVALANCHE_GROWTHRATE ;   R)  @   a   P_STAR%EPAR+CALCULATE_AVALANCHE_GROWTHRATE 8   )  @   a   P_STAR%B+CALCULATE_AVALANCHE_GROWTHRATE 0   Ò)  è       E_CEFF_OVER_E_C+ELECTRIC_FIELDS 8   º*  <      E_CEFF_OVER_E_C%SUM+ELECTRIC_FIELDS=SUM 8   ö*  <      E_CEFF_OVER_E_C%LOG+ELECTRIC_FIELDS=LOG :   2+  =      E_CEFF_OVER_E_C%SQRT+ELECTRIC_FIELDS=SQRT 8   o+  <      E_CEFF_OVER_E_C%ABS+ELECTRIC_FIELDS=ABS 2   «+  ¤   a   E_CEFF_OVER_E_C%Z+ELECTRIC_FIELDS 3   O,  ¤   a   E_CEFF_OVER_E_C%Z0+ELECTRIC_FIELDS 3   ó,  ¤   a   E_CEFF_OVER_E_C%NJ+ELECTRIC_FIELDS 9   -  @   a   E_CEFF_OVER_E_C%NSPECIES+ELECTRIC_FIELDS 2   ×-  @   a   E_CEFF_OVER_E_C%T+ELECTRIC_FIELDS 2   .  @   a   E_CEFF_OVER_E_C%B+ELECTRIC_FIELDS $   W.  _       E_C+ELECTRIC_FIELDS '   ¶.  @   a   E_C%NE+ELECTRIC_FIELDS &   ö.  @   a   E_C%T+ELECTRIC_FIELDS $   6/  _       E_D+ELECTRIC_FIELDS '   /  @   a   E_D%NE+ELECTRIC_FIELDS &   Õ/  @   a   E_D%T+ELECTRIC_FIELDS ,   0  _       NU_EE+COLLISION_FREQUENCIES /   t0  @   a   NU_EE%NE+COLLISION_FREQUENCIES .   ´0  @   a   NU_EE%T+COLLISION_FREQUENCIES .   ô0  ¥       NUBAR_D+COLLISION_FREQUENCIES 6   1  <      NUBAR_D%SUM+COLLISION_FREQUENCIES=SUM 6   Õ1  <      NUBAR_D%LOG+COLLISION_FREQUENCIES=LOG 0   2  ¤   a   NUBAR_D%Z+COLLISION_FREQUENCIES 1   µ2  ¤   a   NUBAR_D%Z0+COLLISION_FREQUENCIES 1   Y3  ¤   a   NUBAR_D%NJ+COLLISION_FREQUENCIES 7   ý3  @   a   NUBAR_D%NSPECIES+COLLISION_FREQUENCIES 0   =4  @   a   NUBAR_D%T+COLLISION_FREQUENCIES 0   }4  @   a   NUBAR_D%P+COLLISION_FREQUENCIES 5   ½4  À       CALC_NUBAR_D01+COLLISION_FREQUENCIES =   }5  <      CALC_NUBAR_D01%SUM+COLLISION_FREQUENCIES=SUM =   ¹5  <      CALC_NUBAR_D01%LOG+COLLISION_FREQUENCIES=LOG 7   õ5  ¤   a   CALC_NUBAR_D01%Z+COLLISION_FREQUENCIES 8   6  ¤   a   CALC_NUBAR_D01%Z0+COLLISION_FREQUENCIES 8   =7  ¤   a   CALC_NUBAR_D01%NJ+COLLISION_FREQUENCIES >   á7  @   a   CALC_NUBAR_D01%NSPECIES+COLLISION_FREQUENCIES 7   !8  @   a   CALC_NUBAR_D01%T+COLLISION_FREQUENCIES >   a8  @   a   CALC_NUBAR_D01%NUBAR_D0+COLLISION_FREQUENCIES >   ¡8  @   a   CALC_NUBAR_D01%NUBAR_D1+COLLISION_FREQUENCIES .   á8  ·       NUBAR_S+COLLISION_FREQUENCIES 6   9  <      NUBAR_S%SUM+COLLISION_FREQUENCIES=SUM 6   Ô9  <      NUBAR_S%LOG+COLLISION_FREQUENCIES=LOG 8   :  =      NUBAR_S%SQRT+COLLISION_FREQUENCIES=SQRT 0   M:  ¤   a   NUBAR_S%Z+COLLISION_FREQUENCIES 1   ñ:  ¤   a   NUBAR_S%Z0+COLLISION_FREQUENCIES 1   ;  ¤   a   NUBAR_S%NJ+COLLISION_FREQUENCIES 7   9<  @   a   NUBAR_S%NSPECIES+COLLISION_FREQUENCIES 0   y<  @   a   NUBAR_S%T+COLLISION_FREQUENCIES 0   ¹<  @   a   NUBAR_S%P+COLLISION_FREQUENCIES 5   ù<  À       CALC_NUBAR_S01+COLLISION_FREQUENCIES =   ¹=  <      CALC_NUBAR_S01%SUM+COLLISION_FREQUENCIES=SUM =   õ=  <      CALC_NUBAR_S01%LOG+COLLISION_FREQUENCIES=LOG 7   1>  ¤   a   CALC_NUBAR_S01%Z+COLLISION_FREQUENCIES 8   Õ>  ¤   a   CALC_NUBAR_S01%Z0+COLLISION_FREQUENCIES 8   y?  ¤   a   CALC_NUBAR_S01%NJ+COLLISION_FREQUENCIES >   @  @   a   CALC_NUBAR_S01%NSPECIES+COLLISION_FREQUENCIES 7   ]@  @   a   CALC_NUBAR_S01%T+COLLISION_FREQUENCIES >   @  @   a   CALC_NUBAR_S01%NUBAR_S0+COLLISION_FREQUENCIES >   Ý@  @   a   CALC_NUBAR_S01%NUBAR_S1+COLLISION_FREQUENCIES /   A  t       LN_LAMBDA_0+COULOMB_LOGARITHMS 7   A  <      LN_LAMBDA_0%LOG+COULOMB_LOGARITHMS=LOG 2   ÍA  @   a   LN_LAMBDA_0%NE+COULOMB_LOGARITHMS 1   B  @   a   LN_LAMBDA_0%T+COULOMB_LOGARITHMS /   MB  t       LN_LAMBDA_C+COULOMB_LOGARITHMS 7   ÁB  <      LN_LAMBDA_C%LOG+COULOMB_LOGARITHMS=LOG 2   ýB  @   a   LN_LAMBDA_C%NE+COULOMB_LOGARITHMS 1   =C  @   a   LN_LAMBDA_C%T+COULOMB_LOGARITHMS 0   }C         LN_LAMBDA_EE+COULOMB_LOGARITHMS 8   D  <      LN_LAMBDA_EE%LOG+COULOMB_LOGARITHMS=LOG :   LD  =      LN_LAMBDA_EE%SQRT+COULOMB_LOGARITHMS=SQRT 3   D  @   a   LN_LAMBDA_EE%NE+COULOMB_LOGARITHMS 2   ÉD  @   a   LN_LAMBDA_EE%T+COULOMB_LOGARITHMS 2   	E  @   a   LN_LAMBDA_EE%P+COULOMB_LOGARITHMS 0   IE         LN_LAMBDA_EI+COULOMB_LOGARITHMS 8   ÜE  <      LN_LAMBDA_EI%LOG+COULOMB_LOGARITHMS=LOG :   F  =      LN_LAMBDA_EI%SQRT+COULOMB_LOGARITHMS=SQRT 3   UF  @   a   LN_LAMBDA_EI%NE+COULOMB_LOGARITHMS 2   F  @   a   LN_LAMBDA_EI%T+COULOMB_LOGARITHMS 2   ÕF  @   a   LN_LAMBDA_EI%P+COULOMB_LOGARITHMS 