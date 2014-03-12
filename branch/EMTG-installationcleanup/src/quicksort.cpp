//quicksort
//algorithm from wikipedia http://en.wikipedia.org/wiki/Quicksort

#include "quicksort.h"
#include <vector>

using namespace std;

namespace EMTG {

//When this is called from an outside function, left and right should be 0 and (size-1) respectively
void quicksort(vector<int>& indices, vector<double>& scores, int left, int right)
{
	 int pivot, leftIdx = left, rightIdx = right;
	 double tempscore;
	 int tempindex;
     if (right - left > 0)
	 {
         pivot = (left + right) / 2;
         while (leftIdx <= pivot && rightIdx >= pivot)
		 {
             while (scores[leftIdx] < scores[pivot] && leftIdx <= pivot)
                 leftIdx = leftIdx + 1;
             while (scores[rightIdx] > scores[pivot] && rightIdx >= pivot)
                 rightIdx = rightIdx - 1;

			 tempscore = scores[rightIdx];
			 scores[rightIdx] = scores[leftIdx];
			 scores[leftIdx] = tempscore;
			 tempindex = indices[rightIdx];
			 indices[rightIdx] = indices[leftIdx];
			 indices[leftIdx] = tempindex;

             leftIdx = leftIdx + 1;
             rightIdx = rightIdx - 1;
             if (leftIdx - 1 == pivot)
                 pivot = rightIdx = rightIdx + 1;
             else if (rightIdx + 1 == pivot)
                 pivot = leftIdx = leftIdx - 1;
		 }

         quicksort(indices, scores, left, pivot - 1);
         quicksort(indices, scores, pivot + 1, right);
	 }
}

}