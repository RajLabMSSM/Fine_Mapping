#ifndef COMBINATION_TABLE_H
#define COMBINATION_TABLE_H

class CombinationTable {
	public:
	CombinationTable(int _max_n, int _max_k) {
		max_n = _max_n;
		max_k = _max_k;
		generateCombinationTable();
	}
	CombinationTable(CombinationTable& other) {
		if (this != &other) {
			max_n = other.max_n;
			max_k = other.max_k;
			generateCombinationTable();
		}
	}
	CombinationTable& operator= (CombinationTable& rhs) {
		if (this != &rhs) {
			max_n = rhs.max_n;
			max_k = rhs.max_k;
			delete [] table;
			generateCombinationTable();
		}
		return *this;
	}
	~CombinationTable() {
		delete [] table;
	}
	static unsigned long long choose(unsigned int n, unsigned int k);
   int generateCombinationTable();
   unsigned long long operator() (unsigned int n, unsigned int k);

    private:
   	int max_n;
	int max_k;
	unsigned long long* table;
};

inline
unsigned long long CombinationTable::operator() (unsigned int n, \
													   unsigned int k) {
	if (n < 1 || n < k) {
		return 0;
	} else {
		return(table[(n - 1) * max_k + (k - 1)]);
	}
}

#endif
