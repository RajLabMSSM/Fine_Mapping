#include "CombinationTable.h"

unsigned long long CombinationTable::choose(unsigned int n, \
							 unsigned int k) {
    if (k > n) {
        return 0;
    }
    unsigned long long r = 1;
    for (unsigned int d = 1; d <= k; ++d) {
        r *= n--;
        r /= d;
    }
    return r;
}

int CombinationTable::generateCombinationTable() {
	table = new unsigned long long[max_k * max_n];
	for (int k = 1; k <= max_k; k++) {
		if (k > max_n) {
		  break;
		}
		for (int n = k; n <= max_n; n++) {
		  table[(n - 1) * max_k + (k - 1)] = choose(n, k);
		}
	}
	return(0);
}


