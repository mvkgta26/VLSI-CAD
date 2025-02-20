#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include "solver.h"

using namespace std;
string input_netlist_textfile = "primary1.txt";  //Holds the name of the input netlist which we want placed on chip surface

class items
{
public:
	bool** netlist;
	bool** gatepadarray;
	float** padarray;
	float** A;
	float** solution;
	float** solutionsorted;
	float** pseudopadarray_for_leftcontainment;
	float** pseudopadarray_for_rightcontainment;
	float* leftgates;
	float* rightgates;
	float* netweight;
	int row;
	int column;
	int numberofpads;
};


void right_containment(items in)
{
	//Create txt file for right containment
	float* rightgates = in.rightgates;
	int halfnumber = in.row-(in.row / 2); //Number of gates / 2
	int numberofpads = in.numberofpads; //Numberofpads on original input netlist
	float** padarray = in.padarray;
	float** pseudopadarray = in.pseudopadarray_for_rightcontainment;


	ofstream out;
	out.open("rightcontainment.txt");


	int num = 1;  //Will be used as gate number on rightside_containment.txt
	for (int i = 0; i < halfnumber; i++)
	{
		ifstream input_netlist;
		input_netlist.open(input_netlist_textfile);
		if (i == 0)   //Execute only on the first iteration
		{
			out << halfnumber << " "; //number of gates (only on the right) printed on out.txt
			out << in.column << endl; //Number of nets printed next on out
		}


		float waster = -1;
		input_netlist >> waster;
		input_netlist >> waster; //Waste the first 2 number on the first row of input_netlist txt
		int n = rightgates[i];
		for (int j = 0; j < n; j++)
		{
			input_netlist >> waster; //Waste the gate number on each row
			int number_of_nets_connected_to_this_gate;
			input_netlist >> number_of_nets_connected_to_this_gate;
			for (int k = 0; k < number_of_nets_connected_to_this_gate; k++)
			{
				input_netlist >> waster;

			}
		}

		//Get the netlist of the actually required gate
		input_netlist >> waster;  //Waste gate number
		out << num << " "; //print gate number on out.txt

		int number_of_nets; //NUmber of nets in required gate
		input_netlist >> number_of_nets;
		out << number_of_nets << " ";
		for (int l = 0; l < number_of_nets; l++)
		{
			int tempnet;
			input_netlist >> tempnet;
			out << tempnet << " ";
		}

		out << endl; //After this gate details are printed, move to next line to print next gate
		input_netlist.close(); //Close so that it can be reopened newly next iteration from first line
		num++;
	}

	//---------------------writing the pads on out.txt

	int number_of_rightcontainmnet_pads = 0;

	for (int i = 0; i < numberofpads; i++)  //Traverse padarray and if x>50 , increment number_of_leftcomtainment_pads
	{
		if (padarray[i][1] >= 50)  //x coorodinate is less than 50
		{
			number_of_rightcontainmnet_pads++;   //Count number of pads on right on surface
		}
	}

	int i = 0;  //For while loop
	while (pseudopadarray[i][0] != -1)   //Count number of rows / pseudo pad gates, in pseudopad array.. (Pseudopad array is filled with -1 in unused cells)
	{
		number_of_rightcontainmnet_pads++;
		i++;
	}

	out << number_of_rightcontainmnet_pads << endl;  //Print numberof left pads in out
	num = 1; //Holds pad number count to be printed in out txt
	for (int i = 0; i < numberofpads; i++)  //Traverse padarray and if x>=50 , print on out.txt
	{
		if (padarray[i][1] >= 50)  //x coorodinate is >= than 50
		{
			out << num << " " << padarray[i][0] << " " << padarray[i][1] << " " << padarray[i][2] << endl;
			num++; //Increment pad count for right containment
		}
	}

	//Write from pseudo pad array on out.txt
	i = 0;
	while (pseudopadarray[i][0] != -1)
	{
		out << num << " " << pseudopadarray[i][1] << " " << pseudopadarray[i][2] << " " << pseudopadarray[i][3] << endl;
		i++;
		num++;  //Increment pad count 
	}

	out.close();


}


void find_gates_in_left_side_connected_to_right_side(items input)
{
	float* leftgates = input.leftgates;
	float* rightgates = input.rightgates;
	float** solutionsorted = input.solutionsorted;
	bool** netlist = input.netlist;
	float** solution = input.solution;
	int numberofgates = input.row;
	int numberofnets = input.column;
	int n1 = numberofgates / 2;  //n1 holds number of elements in leftgates
	int n2 = numberofgates - n1; //n2 holds number of elements in rightgates

	//Creating pseudopadarrayforrightcontainment of rows = n1 and 4 columns: 1 for gatenoindex, 1 for netconnected, 1 for x, 1 for y coordinate
	float** pseudopadarray_for_right_containment = (float**)malloc(n1 * sizeof(float*));
	for (int i = 0; i < n1; i++)
	{
		pseudopadarray_for_right_containment[i] = (float*)malloc(4 * sizeof(float));
	}

	//Initialising all elements in pseudopadarray as -1
	for (int i = 0; i < n1; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			pseudopadarray_for_right_containment[i][j] = -1;

		}
	}

	//-------------------
	int k = 0;

	for (int i = 0; i < n1; i++)	//Traverse through each of elements on leftgates	
	{
		int left_gate_under_comparison = leftgates[i];
		int right_gate_under_comparison;
		bool indicator = 0; //Becomes 1 immediately if the current left gate is found out to have connection to one of right gate by one of nets
		for (int net_number_index = 0; net_number_index < numberofnets; net_number_index++)  //Traverse netlist array through the net_number_index
		{
			for (int j = 0; j < n2; j++)  //Traverse through rightgates elements
			{
				right_gate_under_comparison = rightgates[j];
				//If the present leftgate and rightgate are connected by same net
				if ((netlist[left_gate_under_comparison][net_number_index] == 1) && (netlist[right_gate_under_comparison][net_number_index] == 1))
				{
					pseudopadarray_for_right_containment[k][0] = left_gate_under_comparison;
					pseudopadarray_for_right_containment[k][1] = net_number_index; //Net_number_index connecting both the gates under comparison
					pseudopadarray_for_right_containment[k][2] = solution[left_gate_under_comparison][1]; //Setting the x coordinate of right gate under comparison
					pseudopadarray_for_right_containment[k][3] = solution[left_gate_under_comparison][2]; //Setting the y coordinate of right gate under comparison
					k++;
					indicator = 1;
					break;  //Break the comparison loop of present right gate

				}
			}
			if (indicator == 1)
			{
				break;
			}
		}

	}



	//---------------------------------------------------------Printer------------------------------------------------------------------
	////print pseudopadarray
	//cout << "\nPseudopadarrayforrightcontainment:\n";
	//for (int i = 0; i < n1; i++)
	//{
	//	for (int j = 0; j < 4; j++)
	//	{
	//		cout << pseudopadarray_for_right_containment[i][j] << " ";
	//	}
	//	cout << endl;
	//}


	//-----Propagate pseudo pad gates to centre (Set x=50)----------------------------

	int i = 0;
	while (pseudopadarray_for_right_containment[i][0] != -1)
	{
		pseudopadarray_for_right_containment[i][2] = 50;
		i++;
	}


	input.leftgates = leftgates;
	input.rightgates = rightgates;
	input.pseudopadarray_for_rightcontainment = pseudopadarray_for_right_containment;
	right_containment(input);
}


void left_containment(items in)
{
	//Create txt file for left containment
	float* leftgates = in.leftgates;
	int halfnumber = in.row / 2; //Number of gates / 2
	int numberofpads = in.numberofpads; //Numberofpads on original input netlist
	float** padarray = in.padarray;
	float** pseudopadarray = in.pseudopadarray_for_leftcontainment;


	ofstream out;
	out.open("leftcontainment.txt");


	int num = 1;  //Will be used as gate number on leftside_containment.txt
	for (int i = 0; i < halfnumber; i++)
	{
		ifstream input_netlist;
		input_netlist.open(input_netlist_textfile);
		
		if (i == 0)   //Execute only on the first iteration
		{
			out << halfnumber << " "; //number of gates (only on the left) printed on out.txt
			out << in.column << endl; //Number of nets printed next on out
		}


		float waster = -1;
		input_netlist >> waster;
		input_netlist >> waster; //Waste the first 2 number on the first row of input_netlist txt
		int n = leftgates[i];
		for (int j = 0; j < n; j++)
		{
			input_netlist >> waster; //Waste the gate number on each row
			int number_of_nets_connected_to_this_gate;
			input_netlist >> number_of_nets_connected_to_this_gate;
			for (int k = 0; k < number_of_nets_connected_to_this_gate; k++)
			{
				input_netlist >> waster;

			}
		}

		//Get the netlist of the actually required gate
		input_netlist >> waster;  //Waste gate number
		out << num << " "; //print gate number on out.txt

		int number_of_nets; //NUmber of nets in required gate
		input_netlist >> number_of_nets;
		out << number_of_nets << " ";
		for (int l = 0; l < number_of_nets; l++)
		{
			int tempnet;
			input_netlist >> tempnet;
			out << tempnet << " ";
		}

		out << endl; //After this gate details are printed, move to next line to print next gate
		input_netlist.close(); //Close so that it can be reopened newly next iteration from first line
		num++;
	}

	//---------------------Writing the pads on out txt

	int number_of_leftcontainmnet_pads = 0;

	for (int i = 0; i < numberofpads; i++)  //Traverse padarray and if x<50 , increment number_of_leftcomtainment_pads
	{
		if (padarray[i][1] < 50)  //x coorodinate is less than 50
		{
			number_of_leftcontainmnet_pads++;
		}
	}

	int i = 0;  //For while loop
	while (pseudopadarray[i][0] != -1)   //Count number of rows / pseudo pad gates, in pseudopad array.. (Pseudopad array is filled with -1 in unused cells)
	{
		number_of_leftcontainmnet_pads++;
		i++;
	}

	out << number_of_leftcontainmnet_pads << endl;  //Write numberof left pads in out
	num = 1; //Holds pad number count to be printed in out txt
	for (int i = 0; i < numberofpads; i++)  //Traverse padarray and if x<50 , increment number_of_leftcomtainment_pads
	{
		if (padarray[i][1] < 50)  //x coorodinate is less than 50
		{
			out << num << " " << padarray[i][0] << " " << padarray[i][1] << " " << padarray[i][2] << endl;
			num++; //Increment pad count for left containment
		}
	}

	//Write from pseudo pad array on out.txt
	i = 0;
	while (pseudopadarray[i][0] != -1)
	{
		out << num << " " << pseudopadarray[i][1] << " " << pseudopadarray[i][2] << " " << pseudopadarray[i][3] << endl;
		i++;
		num++;  //Increment pad count 
	}

	out.close();


}


void find_gates_in_right_side_connected_to_left_side(items input)
{
	float* leftgates = input.leftgates;
	float* rightgates = input.rightgates;
	float** solutionsorted = input.solutionsorted;
	bool** netlist = input.netlist;
	float** solution = input.solution;
	int numberofgates = input.row;
	int numberofnets = input.column;
	int n1 = numberofgates / 2;  //n1 holds number of elements in leftgates
	int n2 = numberofgates - n1; //n2 holds number of elements in rightgates

	//Creating pseudopadarrayforleftcontainment of rows = n2 and 4 columns: 1 for gatenoindex, 1 for netconnected, 1 for x, 1 for y coordinate
	float** pseudopadarray_for_left_containment = (float**)malloc(n2 * sizeof(float*));
	for (int i = 0; i < n2; i++)
	{
		pseudopadarray_for_left_containment[i] = (float*)malloc(4 * sizeof(float));
	}

	//Initialising all elements in pseudopadarray as -1
	for (int i = 0; i < n2; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			pseudopadarray_for_left_containment[i][j] = -1;

		}
	}

	//-------------------
	int k = 0;
	

	for (int i = 0; i < n2; i++)	//Traverse through each of elements on rightgates	
	{
		int right_gate_under_comparison = rightgates[i];
		int left_gate_under_comparison;
		bool indicator = 0; //Becomes 1 immediately if the current right gate is found out to have connection to one of left gate by one of nets
		for (int net_number_index = 0; net_number_index < numberofnets; net_number_index++)  //Traverse netlist array through the net_number_index
		{
			for (int j = 0; j < n1; j++)  //Traverse through leftgates elements
			{
				left_gate_under_comparison = leftgates[j];
				//If the present leftgate and rightgate are connected by same net
				if ((netlist[left_gate_under_comparison][net_number_index] == 1) && (netlist[right_gate_under_comparison][net_number_index] == 1))
				{
					pseudopadarray_for_left_containment[k][0] = right_gate_under_comparison;
					pseudopadarray_for_left_containment[k][1] = net_number_index; //Net_number_index connecting both the gates under comparison
					pseudopadarray_for_left_containment[k][2] = solution[right_gate_under_comparison][1]; //Setting the x coordinate of right gate under comparison
					pseudopadarray_for_left_containment[k][3] = solution[right_gate_under_comparison][2]; //Setting the y coordinate of right gate under comparison
					k++;
					indicator = 1;
					break;  //Break the comparison loop of present right gate

				}
			}
			if (indicator == 1)
			{
				break;
			}
		}

	}



	//---------------------------------------------------------Printer------------------------------------------------------------------
	//print pseudopadarray
	cout << "\nPseudopadarrayforledtcontainment:\n";
	for (int i = 0; i < n2; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			cout << pseudopadarray_for_left_containment[i][j] << " ";
		}
		cout << endl;
	}


	//-----Propagate pseudo pad gates to centre (Set x=50)----------------------------

	int i = 0;
	while (pseudopadarray_for_left_containment[i][0] != -1)
	{
		pseudopadarray_for_left_containment[i][2] = 50;
		i++;
	}


	input.leftgates = leftgates;
	input.rightgates = rightgates;
	input.pseudopadarray_for_leftcontainment = pseudopadarray_for_left_containment;
	left_containment(input);
}


void sort_by_x(items in)
{
	float** solution = in.solution;
	int number_of_gates = in.row;
	float* sorted_x = (float*)malloc(number_of_gates * sizeof(float));

	//Creating a dynamic 2D array for solution_sorted
	float** solution_sorted = (float**)malloc(number_of_gates * sizeof(float*));
	for (int i = 0; i < number_of_gates; i++)
	{
		solution_sorted[i] = (float*)malloc(3 * sizeof(float));
	}

	//---------------------
	for (int i = 0; i < number_of_gates; i++)
	{
		sorted_x[i] = solution[i][1];
	}

	//----------------Sorting the solution based on x---------------------------------------
	sort(sorted_x, sorted_x + number_of_gates);



	//---------------------------------------------------------Printer------------------------------------------------------------------
	/*cout << "\n\nSorted_X : " << endl;
	for (int i = 0; i < number_of_gates; i++)
	{
		cout << sorted_x[i] << endl;
	}*/


	int gate_number_index = 0;
	int temp_gate_number_index = 0;  //Take care of equal x repetitions
	for (int i = 0; i < number_of_gates; i++)
	{
		solution_sorted[i][1] = sorted_x[i];
		float temp;

		/*Search the solution array to find the row where x = sorted_x[i] and make solution_sorted[i][2] = y of that corresponding x and
		  solution_sorted[i][0] = gate_number_index of that corresponding x */
		for (int j = 0; j < number_of_gates; j++)
		{
			if (solution[j][1] == sorted_x[i])
			{
				temp = solution[j][2];   //temp holds the y corresponding to current x
				gate_number_index = solution[j][0];   //gate_number_index holds the gate number index of gate whose x coordinate (solution[j][1]) is equal to sorted_x[i]
				if (gate_number_index == temp_gate_number_index) //Take care of repetitions of x
				{
					continue;
				}
				break;
			}
		}
		temp_gate_number_index = gate_number_index;  //need to fix,, more recurrences
		solution_sorted[i][0] = gate_number_index;
		solution_sorted[i][2] = temp;
	}




	//---------------------------------------------------------Printer------------------------------------------------------------------
	////Printing sorted solution
	//cout << "\n\nSorted Solution : " << endl;
	//for (int i = 0; i < number_of_gates; i++)
	//{
	//	cout << solution_sorted[i][0] << " " << solution_sorted[i][1] << " " << solution_sorted[i][2] << endl;
	//}

	int halfnumber = number_of_gates / 2;
	float* leftgates = (float*)malloc(halfnumber * sizeof(float));   //Array of size 9
	float* rightgates = (float*)malloc(halfnumber * sizeof(float));  //Array of size 9

	for (int i = 0; i < halfnumber; i++)
	{
		leftgates[i] = solution_sorted[i][0];
	}

	for (int i = halfnumber; i < number_of_gates; i++)
	{
		rightgates[i - halfnumber] = solution_sorted[i][0];
	}



	//---------------------------------------------------------Printer------------------------------------------------------------------
	/*cout << " \nleftgates: ";
	for (int i = 0; i < halfnumber; i++)
	{
		cout << leftgates[i] << " ";
	}

	cout << " \nrightgates: ";
	for (int i = 0; i < halfnumber; i++)
	{
		cout << rightgates[i] << " ";
	}*/


	//-----------------
	in.leftgates = leftgates;
	in.rightgates = rightgates;
	in.solutionsorted = solution_sorted;
	find_gates_in_right_side_connected_to_left_side(in);
	find_gates_in_left_side_connected_to_right_side(in);

}


void Read_first_solution(items in)
{
	int number_of_gates = in.row;
	int number_of_nets = in.column;

	//Creating a dynamic 2D Array of column size 3 (column 0=gate number index, column1=x, column 2 = y) and row size=number of gates)
	float** solution = (float**)malloc(number_of_gates * sizeof(float*));
	for (int i = 0; i < number_of_gates; i++)
	{
		solution[i] = (float*)malloc(3 * sizeof(float));
	}


	//Reading from solution.txt and storing in solution array
	ifstream sol;
	sol.open("solution.txt");
	int temp, x, y;
	for (int i = 0; i < number_of_gates; i++)
	{
		solution[i][0] = i;
		sol >> temp;
		for (int j = 1; j < 3; j++)
		{
			sol >> solution[i][j];
		}
	}



	//---------------------------------------------------------Printer------------------------------------------------------------------
	////Printing the solution array
	//cout << "\nUnsorted Solution:\n";
	//for (int i = 0; i < number_of_gates; i++)
	//{
	//	for (int j = 0; j < 3; j++)
	//	{
	//		cout << solution[i][j] << " ";
	//	}
	//	cout << endl;
	//}

	in.solution = solution;
	sort_by_x(in);


}


void solve_x_y(items r, ofstream& solution)
{

	int rows = r.row;  //Getting number of rows from r object returned from bmatrices_creator function


	coo_matrix A;
	A.read_coo_matrix("A.txt");
	valarray<double> x(A.n);    //Create an array to implement x matrix
	valarray<double> y(A.n);	//Create an array to implement y matrix
	valarray<double> Bx(A.n);   //Create an array to implement Bx matrix
	valarray<double> By(A.n);   //Create an array to implement By matrix


	ifstream iff("Bx.txt");	//Object iff holds bmatrices_creator.txt file

	//Read from bxmatrix text file and store in Bx array
	for (int i = 0; i < A.n; i++)
	{
		iff >> Bx[i];
	}
	iff.close();


	//Read from bymatrix text file and store in By array
	ifstream iff2("By.txt");
	for (int i = 0; i < A.n; i++)
	{
		iff2 >> By[i];
	}
	iff2.close();




	//---------------------------------------------------------Printer------------------------------------------------------------------
	//cout << "\n" << "Bx Matrix: \n";
	//print_valarray(Bx);  //Print array Bx
	//cout << "\n" << "By Matrix: \n";
	//print_valarray(By);  //Print array By
	A.solve(Bx, x);		 //Solve A*x=Bx
	//cout << "\n\nSolution for x:" << endl;
	//print_valarray(x);    //Print array x
	A.solve(By, y);      //Solve A*y=By  
	//cout << "\n\nSolution for y:" << endl;
	//print_valarray(y);    //Print array y


	//Write solution: "gate number   x     y"   in a file
	//ofstream solution("Solution.txt");
	for (int i = 0; i < rows; i++)
	{
		solution << i + 1 << " " << x[i] << " " << y[i] << " " << endl;
	}
	solution.close();

	//Read_first_solution(r);-----------------------------------------
}


items bmatrices_creator(items r3)
{

	//Getting the passed data from member variables of object r3
	int rows = r3.row;
	int columns = r3.column;
	int numberofpads = r3.numberofpads;
	bool** netlist = r3.netlist;
	float** padarray = r3.padarray;
	bool** gatepadarray = r3.gatepadarray;
	float** A = r3.A;
	float* netweight = r3.netweight;



	//-----------------------------------------------Dynamic Bx matrix--------------------------------------------------------
	float* Bx = (float*)malloc(rows * sizeof(float)); //Create dynamic array Bx of length = number of gates = rows


	//--------------------------------------------Initialising all Bx elements as 0
	for (int i = 0; i < rows; i++)
	{
		Bx[i] = 0;
	}


	//-----------------------------------------------Dynamic By matrix
	float* By = (float*)malloc(rows * sizeof(float)); //Create dynamic array Bx of length = number of gates = rows




	//--------------------------------------------Initialising all By elements as 0----------------------------------------
	for (int i = 0; i < rows; i++)
	{
		By[i] = 0;
	}



	//--------------------------------------------Putting values in Bx matrix------------------------------------------------------
	for (int i = 0; i < rows; i++)			//i holds the gate number
	{
		for (int j = 0; j < numberofpads; j++)  //j iterates through pad numbers
		{
			if (gatepadarray[i][j] == 1)         //If gate i and pad j are connected
			{
				int netconnectedtopad = padarray[j][0] - 1;   //netconnectedtopad holds the net number index of net connected to pad..(Note that padarray[j][0] holds the pad number, not pad number index of the pad_number_index j)
				float temp = netweight[netconnectedtopad];
				Bx[i] = Bx[i] + padarray[j][1] * temp;   //Multiply x coordinate of pad j connected to gate i, and netweight of net connecting these 2 pads.. Add all such multiplications for each pad connection (as j iterates through all of the pads) to gate i 
			}
		}
	}

	//--------------------------------------------Putting values in By matrix------------------------------------
	for (int i = 0; i < rows; i++)			//i holds the gate number
	{
		for (int j = 0; j < numberofpads; j++)  //j iterates through pad numbers
		{
			if (gatepadarray[i][j] == 1)         //If gate i and pad j are connected
			{
				int netconnectedtopad = padarray[j][0] - 1;  //netconnectedtopad holds the net number index of net connected to pad..(Note that padarray[j][0] holds the pad number, not pad number index of the pad_number_index j)
				float temp = netweight[netconnectedtopad];
				By[i] = By[i] + padarray[j][2] * temp;   //Multiply y coordinate of pad j connected to gate i, and netweight of net connecting these 2 pads.. Add all such multiplications for each pad connection (as j iterates through all of the pads) to gate i 
			}
		}
	}


	//----------------------------------------------------Write Bx-------------------------------------------

	//Write elements of Bx to a file named Bx
	ofstream off("Bx.txt");
	for (int i = 0; i < rows; i++)
	{
		off << Bx[i] << endl;
	}
	off.close();

	//Write elements of By to a file named "By.txt"
	ofstream off4("By.txt");
	for (int i = 0; i < rows; i++)
	{
		off4 << By[i] << endl;
	}
	off4.close();



	//Write elements of Matrix to a file named "A.txt"
		//Step 1: Count the number of non zero elements in Matrix (our A matrix)
	int nnzcount = 0;
	for (int i = 0; i < rows; i++)		//Traverse the entire matrix and count the number of non zero elements
	{
		for (int j = 0; j < rows; j++)
		{
			if (A[i][j] != 0)
			{
				nnzcount++;
			}
		}
	}

	ofstream off2("A.txt");
	off2 << rows << " " << nnzcount << endl;   //Write number of gates and number of non zero elements in first line

	//Write non zero elements along with row and column number in A file (One line for one element)

	for (int i = 0; i < rows; i++)		//Traverse the entire matrix and fill in A.txt file as per syntax
	{
		for (int j = 0; j < rows; j++)
		{
			if (A[i][j] != 0)
			{
				off2 << i << " " << j << " " << A[i][j] << endl; //Write: row number i, column number j, non zero element A[i][j] in each line
			}
		}
	}


	return r3;
	//Calling solve function
	//solve_x_y(r3);    //Give r3 as an input object to solve function


}


items A_Matrix_Creator(items r2)
{
	//Getting the passed data from object member variables 
	int rows = r2.row;
	int col = r2.column;
	bool** netlist = r2.netlist;
	float** padarray = r2.padarray;
	int numberofpads = r2.numberofpads;
	bool** gatepadarray = r2.gatepadarray;
	float* netweight = r2.netweight;


	//Creating a dynamic 2D n*n float matrix called Matrix (Rows=n=number of gates).... 
	//Matrix will hold the connectivity matrix initially and later go on to hold the Sparse A matrix

	float** Matrix = (float**)malloc((rows) * sizeof(float*));
	for (int i = 0; i < rows; i++)
	{
		Matrix[i] = (float*)malloc((rows) * sizeof(float));
	}

	//Fill in the Matrix with all 0s
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < rows; j++)
		{
			Matrix[i][j] = 0;
		}
	}




	//Fill appropriate values in the connectivity matrix
	for (int i = 0; i < col; i++)					//i holds the column number or net number index 
	{
		for (int j = 0; j < rows; j++)				//j iterates the row number to catch the first instances of 1 in this column
		{
			if (netlist[j][i] == 1)					//If this row has 1: search for 1s in other rows after this in the same column
			{
				for (int k = j + 1; k < rows; k++)	//k is used to search for second instances of 1 in rows>j
				{
					if (netlist[k][i] == 1)			// If another row has 1, then fill in the 2 cells in Matrix for the correpsonding gate combination with the netweight of netnumber_index
					{
						Matrix[j][k] = netweight[i];
						Matrix[k][j] = netweight[i];
					}
				}
			}
		}
	}


	//-------------------------create sparse matrix from connectivity matrix----------------------------------------------------------

	//Step 1 of connectivity matrix to sparse matrix conversion: make all Matrix elements negative
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < rows; j++)
		{
			Matrix[i][j] = 0 - Matrix[i][j];
		}
	}




	//---------------------------------------------------------Printer------------------------------------------------------------------
	////print modified connectivity matrix
	//cout << "Your connectivity matrix 1st level:\n";
	//for (int i = 0; i < rows; i++)
	//{
	//	for (int j = 0; j < rows; j++)
	//	{
	//		cout << Matrix[i][j] << " ";
	//	}
	//	cout << "\n";
	//}



	//--------------------------------------Setting the Diagonal elements of modified connectivity matrix

	//Step 2 of connectivity matrix to sparse matrix conversion.

	//--------Adding all elements of each row in connectivity matrix and make it the diagonal element of that row in Matrix
	float diagonal = 0;
	/*Iterate through all elements of Matrix row by row and add all elements
	of each row and set as diagonal value for that corresponding row i.e, Matrix[row][row]*/
	for (int i = 0; i < rows; i++)
	{
		diagonal = 0;	//Initialising diagonal=0 after finishing each row's iteration
		for (int j = 0; j < rows; j++)		//Loop to add all elements of row i, store in diagonal
		{
			diagonal = diagonal + (0 - Matrix[i][j]);   //Note: We need to add -Matrix[i][j] because we made Matrix elements negative of what it actually was
		}
		Matrix[i][i] = diagonal;
	}



	//---------------------------------------------------------Printer------------------------------------------------------------------
	////print modified connectivity matrix
	//cout << "\n \nYour connectivity matrix 2nd level: \n";
	//for (int i = 0; i < rows; i++)
	//{
	//	for (int j = 0; j < rows; j++)
	//	{
	//		cout << Matrix[i][j] << " ";
	//	}
	//	cout << "\n";
	//}



	//Step 3 of connectivity matrix to sparse matrix conversion.
	//----------Add to the diagonals at [i][i] the weight of wires coming from pads to the gate i 

	for (int i = 0; i < rows; i++)		//i holds the gate number or Matrix diagonal x=y coordinate
	{
		diagonal = 0;
		for (int j = 0; j < numberofpads; j++)  //j iterates through each pad number_index
		{
			if (gatepadarray[i][j] == 1)
			{
				int tempo = padarray[j][0] - 1;   //Under any problemm replace this in lower line------------------------------------------------------------------------------------------------------------------------------------------------------------
				diagonal = diagonal + netweight[tempo];  //Add to diagonal the weight of net connecting the pad j, that connects with the gate i
			}
		}
		Matrix[i][i] += diagonal;
		/*Adding to the matrix diagonal value at [i][i], the number of pads connected to gate i.
		(Since weight of all nets connecting pads and wires is 1, number of pads connected to each gate is
		same as the sum of weights of all nets from pads that connect the gate)
		*/

	}



	//---------------------------------------------------------Printer------------------------------------------------------------------
	////print modified connectivity matrix final level : the sparse A matrix
	//cout << "\n\nYour connectivity matrix final level ----> sparse A matrix: \n";
	//for (int i = 0; i < rows; i++)
	//{
	//	for (int j = 0; j < rows; j++)
	//	{
	//		cout << Matrix[i][j] << " ";
	//	}
	//	cout << "\n";
	//}



	//Object to be passed on to other functions
	r2.A = Matrix;     // r2.A will hold the sparse A matrix 2D Array pointer
	r2.gatepadarray = gatepadarray;
	r2.netweight = netweight;
	r2 = bmatrices_creator(r2);
	return r2;

}


items Create_Data_Structures_For_Input_Netlist(ifstream& holder)
{
	int arrsize, temp;
	int numberofgates, numberofnets, connectednet;
	/*	ifstream holder;
		holder.open();*/			//Open File
	holder >> numberofgates;			//numberofgates holds the number of gates in input netlist 
	holder >> numberofnets;				//numberofnets holds the number of nets in input netlist 


	

	//Creating 2D array named netlist:
	//A dynamic 2D array of numberofgates rows and numberofnets columns each cell holding a boolean number
	bool** netlist = (bool**)malloc((numberofgates) * sizeof(bool*));
	for (int i = 0; i < numberofgates; i++)
	{
		netlist[i] = (bool*)malloc((numberofnets) * sizeof(bool));
	}

	//Fill in the netlist array with all 0s
	for (int i = 0; i < numberofgates; i++)
	{
		for (int j = 0; j < numberofnets; j++)
		{
			netlist[i][j] = 0;
		}
	}



	//Reading file and storing in netlist array
	for (int i = 0; i < numberofgates; i++)
	{
		holder >> temp;
		holder >> arrsize;
		for (int j = 0; j < arrsize; j++)
		{
			holder >> connectednet;
			netlist[i][connectednet - 1] = 1;
		}
	}



	//---------------------------------------------------------Printer------------------------------------------------------------------
	////Printing the netlist array
	//cout << "\nNetlist:\n";
	//for (int i = 0; i < numberofgates; i++)
	//{
	//	for (int j = 0; j < numberofnets; j++)
	//	{
	//		cout << netlist[i][j] << " ";
	//	}
	//	cout << endl;
	//}




	//Reading pad data and storing into a new array called padarray
	int numberofpads = 0;
	holder >> numberofpads;  //Read number of pads

	//Each column is a net and each row is a pad. If connected 1 otherwise 0.
	//Creating dynamic 2D Array named padarray for storing pad data
	float** padarray = (float**)malloc(numberofpads * sizeof(float*));
	for (int i = 0; i < numberofpads; i++)
	{
		padarray[i] = (float*)malloc(3 * sizeof(int));
	}




	//Storing from file into padarray
	for (int i = 0; i < numberofpads; i++)
	{
		holder >> temp;
		for (int j = 0; j < 3; j++)
		{
			holder >> padarray[i][j];
		}
	}

	holder.close();					//Close file


//---------------------------------------------------------Printer------------------------------------------------------------------
	////Printing padarray elements
	//cout << "\n\nPadarray elements: \n";
	//for (int i = 0; i < numberofpads; i++)
	//{
	//	for (int j = 0; j < 3; j++)
	//	{
	//		cout << padarray[i][j] << " ";
	//	}
	//	cout << "\n";
	//}



	//---------------------------------------------------Gatepadarray-----------------------------------------------------------------
	//Gate and pad connection array.....Each row represents a gate and each column represents a pad. Each element is 1 if that gate and the corresponding pad are connected
	//int numberofpads = r2.numberofpads;

	//Creating a dynamic 2D gate pad array called "gatepadarray"
	bool** gatepadarray = (bool**)malloc(numberofgates * sizeof(bool*));
	for (int i = 0; i < numberofgates; i++)
	{
		gatepadarray[i] = (bool*)malloc(numberofpads * sizeof(bool));
	}


	//Make all elements in gatepadarray default as 0 
	for (int i = 0; i < numberofgates; i++)
	{
		for (int j = 0; j < numberofpads; j++)
		{
			gatepadarray[i][j] = 0;
		}

	}


	//Fill the gatepadarray with appropriate values
	int netconnected;
	//This nested loop checks if each gate is connected to same net as each pad. If yes: That means this gate and this pad are connected: Hence the corresponding element in gatepadarray is made 1
	for (int i = 0; i < numberofgates; i++)	//At any itertation: i holds the gate number(row number) and j holds the pad number(column number)
	{
		for (int j = 0; j < numberofpads; j++)
		{
			netconnected = padarray[j][0];  //Is the net that is connected by pad j
			if (netlist[i][netconnected - 1] == 1)	//Netlist[i][netconnected-1] is 1 if gate i is connected to the netconnected net 
			{
				gatepadarray[i][j] = 1;
			}
		}
	}




	//---------------------------------------------------------Printer------------------------------------------------------------------
	////Print gatepad array
	//cout << "\n\nGatepadarray: \n";
	//for (int i = 0; i < numberofgates; i++)
	//{
	//	for (int j = 0; j < numberofpads; j++)
	//	{
	//		cout << gatepadarray[i][j] << " ";
	//	}
	//	cout << "\n";
	//}
	//cout << "\n\n";



	//------------------------------------------------Netweight calculator 
	//( All the weights are actually 1 in input netlist, but we assume new weights of nets according to the number of points they connect to)

	//Create a dynamic array of length = number of nets .  netweight[i] holds the weight of the net number_index i 
	float* netweight = (float*)malloc(numberofnets * sizeof(float));

	for (int i = 0; i < numberofnets; i++)   //Initialise all elements of netweight array as 0
	{
		netweight[i] = 0;
	}

	for (int i = 0; i < numberofnets; i++)			  //i holds the column number or net number index
	{
		int netsize = 0;			//Reset netsize for each new column
		float temp = 1.0;
		for (int j = 0; j < numberofgates; j++)			 //j traverses the rows 
		{
			if (netlist[j][i] == 1)			// If row j has 1, then netsize++
			{
				netsize++;                 //Netsize is number of points (or gates) that the net connects.. Note: We did not consider the nets' connection to pads.
			}
		}


		//Traverse column 0 of padarray and count number of pads connected to net i and increment netsize
		for (int j = 0; j < numberofpads; j++)
		{
			if (padarray[j][0] == i + 1)  // i+1 because we are converting net number index to net number... This is because column 1 of padarray holds net number and not net index.
			{
				netsize++;
			}
		}

		//Netsize now holds the total number of gates and pads that net i connects

		temp = 1.0 / (netsize - 1);    //temp = netweight = 1/(k-1) for a k point net
		netweight[i] = temp;    //Store netweight of net number_index i at netweight[i]

	}




	//-----------------------------------------------------------------------------------------------------------------------
	//ret object holds as member variables the data parameters that are to be  passed on to other functions
	items ret;
	ret.column = numberofnets;
	ret.row = numberofgates;
	ret.netlist = netlist;
	ret.padarray = padarray; //Storing padarray pointer in ret object
	ret.numberofpads = numberofpads; //Storing numberofpads in ret object
	ret.gatepadarray = gatepadarray;   //Storing gatepadarray pointer in r2 object
	ret.netweight = netweight;


	ret = A_Matrix_Creator(ret);
	return ret;
}


void merge_left_and_right_solutions(int numberofgates)
{
	ifstream inp1;
	inp1.open("leftsolution.txt");
	ifstream inp2;
	inp2.open("rightsolution.txt");
	ofstream out;
	out.open("QP3Solution.txt");
	int num = 1;
	float temp;
	for (int i=0; i<(numberofgates/2); i++) //Iterate as many times as number of gates in left
	{
		inp1 >> temp; //waste the gate number on each line available on left and right solutions, make our own count
		out << num;  //Print our own gate count
		for (int i = 0; i < 2; i++) //Read twice from each line and write
		{
			inp1 >> temp;
			out << " " << temp << " ";
			
		}
		num++; //Increment gate number
		out << endl; //Move to next line in write file
	}

	inp1.close();

	for (int i = 0; i < (numberofgates-(numberofgates / 2)); i++)    //Iterate as many times as number of gates in right 
	{
		inp2 >> temp; //waste the gate number on each line available on left and right solutions, make our own count
		out << num;  //Print our own gate count
		for (int i = 0; i < 2; i++) //Read twice from each line and write
		{
			inp2 >> temp;
			out << " " << temp << " ";

		}
		num++; //Increment gate number
		out << endl; //Move to next line in write file

	}

	inp2.close();
	out.close();

}

void main()
{
	ifstream holder;
	(holder).open(input_netlist_textfile);
	items items_for_QP0 =	Create_Data_Structures_For_Input_Netlist(holder);
	int number_of_gates = items_for_QP0.row;  //Number of gates in the main input netlist
	ofstream solutionfile;
	solutionfile.open("solution.txt");
	solve_x_y(items_for_QP0, solutionfile);
	solutionfile.close();
	Read_first_solution(items_for_QP0);  //Read QP0 and create leftcontainment.txt and rightcontainment.txt file


	//Passing the left containment solution to the Placer recursively
	ifstream input_for_leftcontainment;
	input_for_leftcontainment.open("leftcontainment.txt");
	items left_items = Create_Data_Structures_For_Input_Netlist(input_for_leftcontainment);
	ofstream left_solution;
	left_solution.open("leftsolution.txt");
	solve_x_y(left_items, left_solution);


	//Passing the right containment solution to the Placer recursively
	ifstream input_for_rightcontainment;
	input_for_rightcontainment.open("rightcontainment.txt");
	items right_items = Create_Data_Structures_For_Input_Netlist(input_for_rightcontainment);
	ofstream right_solution;
	right_solution.open("rightsolution.txt");
	solve_x_y(right_items, right_solution);



	merge_left_and_right_solutions(number_of_gates);


}
