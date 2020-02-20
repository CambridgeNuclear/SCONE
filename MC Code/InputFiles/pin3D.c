type eigenPhysicsPackage; 

//pop      100000; 
pop      20000;
active   500; 
inactive 300; 
//seed     -6574747;
XSdata   ce2;

collisionOperator { neutronCE {type neutronCEstd;} 
                  } 

transportOperator { type transportOperatorST;
                    // type transportOperatorDT;     
                  } 

inactiveTally { //clerk1 { type keffActiveClerk; display yes; }
                SEClerkIn {type	shannonEntropyClerk;
			   map { type spaceMap; grid lin; axis z; min -183.0; max 183.0; N 20 ;}
			   cycles 300;}
		comClerkIn {type centreOfMassClerk;
			    cycles 300;}
              } 

activeTally  {
	       SEClerkA {type shannonEntropyClerk;
		                                  map {type spaceMap; grid lin; axis z; min -183.0; max 183.0; N 20;}
						  cycles 500;}
	       comClerkA {type centreOfMassClerk;
		          cycles 500;}
	      }

geometry { 
    type basicCellCSG;
    boundary (1 1 1 1 0 0);

    surfaces
    {
      zcyl        { id 2; type zCylinder      ; origin ( 0.0  0.0  0.0); radius 0.4095;}
      box         { id 1; type box            ; origin ( 0.0  0.0  0.0); halfwidth (0.63 0.63 183.0);}
    }


    cells
    {

      fuel {id 1; surfaces (-2); filltype mat; mat uo2;  }
      mod  {id 2; surfaces (2) ; filltype mat; mat water; }
  
  
      out     { id 3; surfaces (1 ); filltype outside;         }
      inside  { id 4; surfaces (-1); filltype uni; universe 20;}
    }

    universes
    {

      root
      {
	  id 1;
	  type cellUniverse; 
	  origin (0.0 0.0 0.0);
	  cells ( 4 3); 
      }

      uni20 {id 20; type cellUniverse; origin (0.0 0.0 0.0); cells (1 2); }     
    }   
} 


nuclearData {

  handles { 
      ce2 byNucMT; 
  }
  
  
materials { 
    aceLibrary /home/pmc55/myACE/JEF311.aceXS; 

    water { 
      temp      273; 
      1001.03c    5.028E-02;
      8016.03c    2.505E-02;
      5010.03c    2.0E-005;
      //xsFile ./InputFiles/CAS8_12_1;  
      xsFile ./InputFiles/WIMS69HS_12_1;
      
    } 

    clad {
      temp 273;
      40090.03c  2.224658E-02;
      40091.03c  0.485144E-02;
      40092.03c  0.741553E-02;
      40094.03c  0.751498E-02;
      40096.03c  0.121070E-02;     
      xsFile ./InputFiles/WIMS69HS_12_1;
    }
    
    
    uo2 {  
      temp       273; 
      8016.03c   4.571E-02; 
      92238.03c  2.180E-02; 
      92235.03c  1.042E-03;
      //xsFile ./InputFiles/CAS8_11_1;
      xsFile ./InputFiles/WIMS69HS_11_1;      
    }	 

} 
  
}
  
  
  
