import java.util.*;
import java.io.FileWriter;
public class ComparingSortingAlgorithms {
    
   static void insertionSort(int arr[]){
		
		int n = arr.length;
        for (int i = 1; i < n; ++i) {
            int key = arr[i];
            int j = i - 1;
 
            /* Move elements of arr[0..i-1], that are
               greater than key, to one position ahead
               of their current position */
            while (j >= 0 && arr[j] > key) {
               
            	arr[j + 1] = arr[j];
                j = j - 1;
            }
            arr[j + 1] = key;
        }
	}  
   
   static void binaryInsertionSort(int nums[])
   {
       int count = 0;
       for (int i = 1; i < nums.length; ++i) {
           int key = nums[i];
           int insertedPosition = findPosition(nums, 0, i - 1, key);

           for (int j = i - 1; j >= insertedPosition; --j) {
               nums[j + 1] = nums[j];
           }

           nums[insertedPosition] = key;
       }
		
   }


	static int findPosition(int[] nums, int start, int end, int key) {
	    while (start <= end) {
	        int mid = start + (end - start) / 2;

	        if (key < nums[mid]) {
	            end = mid - 1;
	        } else {
	            start = mid + 1;
	        }
	    }

	    return start;
	}
	
	
	
	static void merge(int Arr[], int start, int mid, int end) {

		// create a temp array
		int temp[] = new int[end - start + 1];

		// crawlers for both intervals and for temp
		int i = start, j = mid+1, k = 0;

		// traverse both arrays and in each iteration add smaller of both elements in temp 
		while(i <= mid && j <= end) {
			if(Arr[i] <= Arr[j]) {
				temp[k] = Arr[i];
				k += 1; i += 1;
			}
			else {
				temp[k] = Arr[j];
				k += 1; j += 1;
			}
		}

		// add elements left in the first interval 
		while(i <= mid) {
			temp[k] = Arr[i];
			k += 1; i += 1;
		}

		// add elements left in the second interval 
		while(j <= end) {
			temp[k] = Arr[j];
			k += 1; j += 1;
		}

		// copy temp to original interval
		for(i = start; i <= end; i += 1) {
			Arr[i] = temp[i - start];
		}
	}

	static void mergeSort(int Arr[], int start, int end) {

		if(start < end) {
			int mid = (start + end) / 2;
			mergeSort(Arr, start, mid);
			mergeSort(Arr, mid+1, end);
			merge(Arr, start, mid, end);
		}
	}
	
	
	//Function that implements the quick sort algorithm by selecting first element as pivot.
    static void quick_sort_first(int arr[],int low,int high) {
	 
	    int pivot=arr[low];   // Pivot is selected as first element.
	    //Then, Hoare's partitioning algorithm is executed.
	    int i=low-1;
	    int j=high+1;
	    int temp;
	    while(i < j) { 
		
		    do {
			    i++;
		    }while(arr[i] < pivot);
		
		    do {
			    j--;
		    }while(arr[j ] > pivot);
	    
		    if(i<j) {
			    temp=arr[j];  arr[j]=arr[i];  arr[i]=temp; }
	        }
	
       temp=arr[low];  arr[low]=arr[j];  arr[j]=temp;
       
       //After partitioning, recursive statements are executed to sort sub-arrays.
	   if(low<high) {	
		   quick_sort_first(arr,low,j);
		   quick_sort_first(arr,j+1,high);    		
	   }
	
   }
    
    //Function that implements the quick sort algorithm with median-of-three pivot selection.
    static void quick_sort_median_of_three(int arr[],int low,int high) {
  	    
    	//First median of three is found and assigned to pivot.
    	int middle = low + (high-low) / 2;
    	int pivot;
        
        if ((arr[low] > arr[middle]) ^ (arr[low] > arr[high])) 
            pivot = arr[low];
        else if ((arr[middle] < arr[low]) ^ ( arr[middle] < arr[high])) 
        	pivot = arr[middle];
        else
        	pivot = arr[high];
    	
        //Then, Hoare's partitioning is implemented.
        int i=low-1;
    	int j=high+1;
    	int temp;
    	while(i < j) {
    		
    		do {
    			i++;
    		}while(arr[i] < pivot);
    		
    		do {
    			j--;
    		}while(arr[j ] > pivot);
    	    
    		if(i<j) {
    			 temp=arr[j];  arr[j]=arr[i];  arr[i]=temp; }
    	}
    	
        temp=arr[low];  arr[low]=arr[j];  arr[j]=temp;
        
        //After partition, recursive statements are executed to sort sub-arrays.
    	if(low<high) {	
    		quick_sort_median_of_three(arr,low,j);
    		quick_sort_median_of_three(arr,j+1,high);    		
    	}
    	
    }
    
    
    
    public static void heapSort(int arr[])
    {
        int n = arr.length;
 
        // Build heap (rearrange array)
        for (int i = n / 2 - 1; i >= 0; i--)
            heapify(arr, n, i);
 
        // One by one extract an element from heap
        for (int i = n - 1; i > 0; i--) {
            // Move current root to end
            int temp = arr[0];
            arr[0] = arr[i];
            arr[i] = temp;
 
            // call max heapify on the reduced heap
            heapify(arr, i, 0);
        }
    }
    
    static void heapify(int arr[], int n, int i)
    {
        int largest = i; // Initialize largest as root
        int l = 2 * i + 1; // left = 2*i + 1
        int r = 2 * i + 2; // right = 2*i + 2
 
        // If left child is larger than root
        if (l < n && arr[l] > arr[largest])
            largest = l;
 
        // If right child is larger than largest so far
        if (r < n && arr[r] > arr[largest])
            largest = r;
 
        // If largest is not root
        if (largest != i) {
            int swap = arr[i];
            arr[i] = arr[largest];
            arr[largest] = swap;
 
            // Recursively heapify the affected sub-tree
            heapify(arr, n, largest);
        }
    }
    
    
    
    
    public static void countSort(int array[], int size) {
        int[] output = new int[size + 1];

        // Find the largest element of the array
        int max = array[0];
        for (int i = 1; i < size; i++) {
          if (array[i] > max)
            max = array[i];
        }
        int[] count = new int[max + 1];

        // Initialize count array with all zeros.
        for (int i = 0; i < max; ++i) {
          count[i] = 0;
        }

        // Store the count of each element
        for (int i = 0; i < size; i++) {
          count[array[i]]++;
        }

        // Store the cumulative count of each array
        for (int i = 1; i <= max; i++) {
          count[i] += count[i - 1];
        }

        // Find the index of each element of the original array in count array, and
        // place the elements in output array
        for (int i = size - 1; i >= 0; i--) {
          output[count[array[i]] - 1] = array[i];
          count[array[i]]--;
        }

        // Copy the sorted elements into original array
        for (int i = 0; i < size; i++) {
          array[i] = output[i];
        }
      }
    
    //Function that takes a sorted array and gives its reverse sorted version. It will be used to generate reversed sorted array.
      static int[] reverse(int a[], int n) {
    
        int[] b = new int[n];
        int j = n;
        for (int i = 0; i < n; i++) {
            b[j - 1] = a[i];
            j = j - 1;
         }
	  return b;
    
      }

      
    //Main Function
    public static void main(String[] args ) throws Exception {
		
    	Random rand = new Random();
    	int element; //Holds elements of each array.
    	int size=1000;  // Holds size of the arrays.
    	long tStart;  //Holds starting time for an algorithm.
    	long tFinish;  //Holds finishing time for an algorithm
		long tReal;    // Holds start-finish, real time consumed by an algorithm
		long average = 0;  // Holds average time.
    	FileWriter writer = new FileWriter("inputs.txt");
    	
    	/*Using FileWriter class, we create a "inputs.txt" file including the inputs. We first generate input arrays then write
    	* them to input file. In each execution of code, different input files are created but all have same types of inputs.
    	* We submitted one of inputs that generated and we prepared the report according to that input file. But when the code is 
    	* executed again different input file also will be created.*/
   
    	//Three types of inputs will be generated, first is random.
    	//The brief information will be given on console for all types of inputs and all sizes from 1k to 10k. 	
    	//First random arrays.
    	System.out.println("--------------Random elements--------------" );
    	//Following loop generates random arrays for changing sizes, its elements are between 0-10000.
    	while(size<=10000) {    //Size from 1k to 10k increasing by 1k.
    	    
    		int arr [] = new int [size];   	    
    		for(int i = 0; i<size; i++) {
    		    element=rand.nextInt(10001);
    		    arr[i]=element;
    	    }
            
            //This random array is written to inputs.txt.
    		writer.write("Random elements of size " + size + ": \n" );  		
    		for(int i = 0; i<size; i++) {
        	    writer.write(arr[i] + " " );
        	    if(i>0 && i%50==0)
        	    	writer.write("\n");
            }			
    		writer.write("\n\n\n");
    		   		
    		System.out.println("size: " + size );   // Size information.
    		
    		//This and other print functions inform to user on console.
    		System.out.print("Insertion Sort - Time Complexity:O(n^2) ---- ");		
    		//Insertion sort is executed five times and average is taken.
    		for(int x = 0 ; x<5 ; x++) {	
      		  tStart=System.nanoTime();
      		  insertionSort(arr);
              tFinish=System.nanoTime();
      		  tReal = tFinish-tStart;
      		  average += tReal; 
      		}
    		average = average/5;   //average is taken as nanoseconds.
    		average = average/1000; //Nanoseconds is transformed to microseconds.
    		System.out.println("Time spent on average is: " + average + " micro-seconds.");
    		average = 0;    // Average is zero again to perform next algorithm.
    		
    		//Same principle is maintained.
    		System.out.print("Binary Insertion Sort - Time Complexity:O(n^2) ---- ");		//Time complexity for such a input.
    		for(int x = 0 ; x<5 ; x++) {	
      		  tStart=System.nanoTime();
      		  binaryInsertionSort(arr);
              tFinish=System.nanoTime();
      		  tReal = tFinish-tStart;
      		  average += tReal; 
      		}
    		average = average/5;
    		average = average/1000;
    		System.out.println("Time spent on average is: " + average + " micro-seconds." );
    		average = 0;
    		
    		//Merge Sort.
    		//These complexities also can be given as big teta notation but java doesn't include that character that's why, we used big oh notation.
    		System.out.print("Merge-Sort - Time Complexity:O(nlogn) ---- "); 	
    		for(int x = 0 ; x<5 ; x++) {	
      		  tStart=System.nanoTime();
      		  mergeSort(arr , 0 , arr.length-1);
              tFinish=System.nanoTime();
      		  tReal = tFinish-tStart;
      		  average += tReal; 
      		}
    		average = average/5;
    		average = average/1000;
    		System.out.println("Time spent on average is: " + average + " micro-seconds." );
    		average = 0;
    		
    		//Quick sort with first pivot selection
    		System.out.print("Quick Sort with first pivot selection - Time Complexity:O(nlogn) ---- ");		
    		for(int x = 0 ; x<5 ; x++) {	
      		  tStart=System.nanoTime();
      		  quick_sort_first(arr, 0 , arr.length-1);
              tFinish=System.nanoTime();
      		  tReal = tFinish-tStart;
      		  average += tReal; 
      		}
    		average = average/5;
    		average = average/1000;
    		System.out.println("Time spent on average is: " + average + " micro-seconds." );
    		average = 0;
    		
    		//Quick sort with median of three pivot selection.
    		System.out.print("Quick-Sort with median-of-three pivot selection - Time Complexity:O(nlogn) ---- ");		
    		for(int x = 0 ; x<5 ; x++) {	
      		  tStart=System.nanoTime();
      		  quick_sort_median_of_three(arr, 0 , arr.length-1);
              tFinish=System.nanoTime();
      		  tReal = tFinish-tStart;
      		  average += tReal; 
      		}
    		average = average/5;
    		average = average/1000;
    		System.out.println("Time spent on average is: " + average + " micro-seconds." );
    		average = 0;
    		
    		//Heap sort
    		System.out.print("Heap-Sort - Time Complexity:O(nlogn) ---- ");		
    		for(int x = 0 ; x<5 ; x++) {	
      		  tStart=System.nanoTime();
      		  heapSort(arr);
              tFinish=System.nanoTime();
      		  tReal = tFinish-tStart;
      		  average += tReal; 
      		}
    		average = average/5;
    		average = average/1000;
    		System.out.println("Time spent on average is: " + average + " micro-seconds.");
    		average = 0;
    	
    		//Count Sort
    		System.out.print("Count-Sort - Time Complexity:O(n) ---- ");		
    		for(int x = 0 ; x<5 ; x++) {	
      		  tStart=System.nanoTime();
      		  countSort(arr,size);
              tFinish=System.nanoTime();
      		  tReal = tFinish-tStart;
      		  average += tReal; 
      		}
    		average = average/5;
    		average = average/1000;
    		System.out.println("Time spent on average is: " + average + " micro-seconds.");
    		average = 0;
    	    
           System.out.println();    		
    	   size +=1000;           // Size is increased by 1k in each iteration.
    		
    	}
    	
    	size = 1000;   // Size is again 1000 because new type of input will be generated and evaluated now.
    	
    	//Sorted arrays from 1k to 10k.
    	System.out.println("--------------Sorted elements--------------" );
    	while(size <= 10000) {
    		
    		//First random arrays are generated.
    		int arr [] = new int [size];   	    
    		for(int i = 0; i<size; i++) {
    		    element=rand.nextInt(10001);
    		    arr[i]=element;
    	    }
    		 		
    		writer.write("Sorted elements of size " + size + ": \n" );    
    		// They will be sorted and written to "inputs" file.	 
    		Arrays.sort(arr);   		 
    		for(int i = 0; i<size; i++) {
        	    writer.write(arr[i] + " " );
        	    if(i>0 && i%50==0)
        	    	writer.write("\n");
            }    		
    		writer.write("\n\n\n");
    		  		
    		System.out.println("size: " + size);  //Size information.
    		
    		//All types will be generated and user will be informed again with same procedure.
    		System.out.print("Insertion Sort - Time Complexity:O(n) ---- ");		
    		for(int x = 0 ; x<5 ; x++) {	
      		  tStart=System.nanoTime();
      		  insertionSort(arr);
              tFinish=System.nanoTime();
      		  tReal = tFinish-tStart;
      		  average += tReal; 
      		}
    		average = average/5;
    		average = average/1000;
    		System.out.println("Time spent on average is: " + average + " micro-seconds.");
    		average = 0;
    		
    		System.out.print("Binary Insertion Sort - Time Complexity:O(nlogn) ---- ");		
    		for(int x = 0 ; x<5 ; x++) {	
      		  tStart=System.nanoTime();
      		  binaryInsertionSort(arr);
              tFinish=System.nanoTime();
      		  tReal = tFinish-tStart;
      		  average += tReal; 
      		}
    		average = average/5;
    		average = average/1000;
    		System.out.println("Time spent on average is: " + average + " micro-seconds." );
    		average = 0;
    		
    		System.out.print("Merge-Sort - Time Complexity:O(nlogn) ---- ");		
    		for(int x = 0 ; x<5 ; x++) {	
      		  tStart=System.nanoTime();
      		  mergeSort(arr , 0 , arr.length-1);
              tFinish=System.nanoTime();
      		  tReal = tFinish-tStart;
      		  average += tReal; 
      		}
    		average = average/5;
    		average = average/1000;
    		System.out.println("Time spent on average is: " + average + " micro-seconds." );
    		average = 0;
    		
    		System.out.print("Quick Sort with first pivot selection - Time Complexity:O(n^2) ---- ");		
    		for(int x = 0 ; x<5 ; x++) {	
      		  tStart=System.nanoTime();
      		  quick_sort_first(arr, 0 , arr.length-1);
              tFinish=System.nanoTime();
      		  tReal = tFinish-tStart;
      		  average += tReal; 
      		}
    		average = average/5;
    		average = average/1000;
    		System.out.println("Time spent on average is: " + average + " micro-seconds." );
    		average = 0;
    		
    		
    		System.out.print("Quick-Sort with median-of-three pivot selection - Time Complexity:O(nlogn) ---- ");		
    		for(int x = 0 ; x<5 ; x++) {	
      		  tStart=System.nanoTime();
      		  quick_sort_median_of_three(arr, 0 , arr.length-1);
              tFinish=System.nanoTime();
      		  tReal = tFinish-tStart;
      		  average += tReal; 
      		}
    		average = average/5;
    		average = average/1000;
    		System.out.println("Time spent on average is: " + average + " micro-seconds." );
    		average = 0;
    		
    		System.out.print("Heap-Sort - Time Complexity:O(nlogn) ---- ");		
    		for(int x = 0 ; x<5 ; x++) {	
      		  tStart=System.nanoTime();
      		  heapSort(arr);
              tFinish=System.nanoTime();
      		  tReal = tFinish-tStart;
      		  average += tReal; 
      		}
    		average = average/5;
    		average = average/1000;
    		System.out.println("Time spent on average is: " + average + " micro-seconds.");
    		average = 0;
    		
    		System.out.print("Count-Sort - Time Complexity:O(n) ---- ");		
    		for(int x = 0 ; x<5 ; x++) {	
      		  tStart=System.nanoTime();
      		  countSort(arr,size);
              tFinish=System.nanoTime();
      		  tReal = tFinish-tStart;
      		  average += tReal; 
      		}
    		
    		average = average/5;
    		average = average/1000;
    		System.out.println("Time spent on average is: " + average + " micro-seconds.");
    		average = 0;
    	    
    	    System.out.println();
    		size+=1000;
    	    
    	}
    	
    	size=1000;
    	
    	
    	//Reverse Sorted Elements
    	System.out.println("--------------Reverse sorted elements--------------" );
    	while(size <= 10000) {
    		
    		//Random arrays is generated again for all sizes.
    		int arr [] = new int [size];   	    
    		for(int i = 0; i<size; i++) {
    		    element=rand.nextInt(10001);
    		    arr[i]=element;
    	    }    		
    		
    		writer.write("Reverse sorted elements of size " + size + ": \n" );   		
    		//Arrays will be sorted first and just after reverse sorted then written to "inputs" file.
    		Arrays.sort(arr);
    		arr = reverse(arr,size);   		
    		for(int i = 0; i<size; i++) {
        	    writer.write(arr[i] + " " );
        	    if(i>0 && i%50==0)
        	    	writer.write("\n");
            }	 
    		writer.write("\n\n\n");  		
    		
    		System.out.println("size: " + size);
    		
    		//All input types for changing size.
    		System.out.print("Insertion Sort - Time Complexity:O(n^2) ---- ");		
    		for(int x = 0 ; x<5 ; x++) {	
      		  tStart=System.nanoTime();
      		  insertionSort(arr);
              tFinish=System.nanoTime();
      		  tReal = tFinish-tStart;
      		  average += tReal; 
      		}
    		average = average/5;
    		average = average/1000;
    		System.out.println("Time spent on average is: " + average + " micro-seconds.");
    		average = 0;
    		
    		System.out.print("Binary Insertion Sort - Time Complexity:O(n^2) ---- ");		
    		for(int x = 0 ; x<5 ; x++) {	
      		  tStart=System.nanoTime();
      		  binaryInsertionSort(arr);
              tFinish=System.nanoTime();
      		  tReal = tFinish-tStart;
      		  average += tReal; 
      		}
    		average = average/5;
    		average = average/1000;
    		System.out.println("Time spent on average is: " + average + " micro-seconds." );
    		average = 0;
    		
    		System.out.print("Merge-Sort - Time Complexity:O(nlogn) ---- ");		
    		for(int x = 0 ; x<5 ; x++) {	
      		  tStart=System.nanoTime();
      		  mergeSort(arr , 0 , arr.length-1);
              tFinish=System.nanoTime();
      		  tReal = tFinish-tStart;
      		  average += tReal; 
      		}
    		average = average/5;
    		average = average/1000;
    		System.out.println("Time spent on average is: " + average + " micro-seconds." );
    		average = 0;
    		
    		System.out.print("Quick Sort with first pivot selection - Time Complexity:O(n^2) ---- ");		
    		for(int x = 0 ; x<5 ; x++) {	
      		  tStart=System.nanoTime();
      		  quick_sort_first(arr, 0 , arr.length-1);
              tFinish=System.nanoTime();
      		  tReal = tFinish-tStart;
      		  average += tReal; 
      		}
    		average = average/5;
    		average = average/1000;
    		System.out.println("Time spent on average is: " + average + " micro-seconds." );
    		average = 0;
    		
    		
    		System.out.print("Quick-Sort with median-of-three pivot selection - Time Complexity:O(nlogn) ---- ");		
    		for(int x = 0 ; x<5 ; x++) {	
      		  tStart=System.nanoTime();
      		  quick_sort_median_of_three(arr, 0 , arr.length-1);
              tFinish=System.nanoTime();
      		  tReal = tFinish-tStart;
      		  average += tReal; 
      		}
    		average = average/5;
    		average = average/1000;
    		System.out.println("Time spent on average is: " + average + " micro-seconds." );
    		average = 0;
    		
    		System.out.print("Heap-Sort - Time Complexity:O(nlogn) ---- ");		
    		for(int x = 0 ; x<5 ; x++) {	
      		  tStart=System.nanoTime();
      		  heapSort(arr);
              tFinish=System.nanoTime();
      		  tReal = tFinish-tStart;
      		  average += tReal; 
      		}
    		average = average/5;
    		average = average/1000;
    		System.out.println("Time spent on average is: " + average + " micro-seconds.");
    		average = 0;
    		
    		System.out.print("Count-Sort - Time Complexity-O(n) ---- ");		
    		for(int x = 0 ; x<5 ; x++) {	
      		  tStart=System.nanoTime();
      		  countSort(arr,size);
              tFinish=System.nanoTime();
      		  tReal = tFinish-tStart;
      		  average += tReal; 
      		}
    		average = average/5;
    		average = average/1000;
    		System.out.println("Time spent on average is: " + average + " micro-seconds.");
    		average = 0;
    		
    		System.out.println();  		
            size+=1000;
    	}	
    	
    	writer.close();
    	  	
    
	}    
} 
