#pragma once

#include <stdint.h>
#include <intrin.h>

template<typename T>
class Heap
{
public:

	Heap() :m_innerHeap(0), m_allocated(0), m_size(0) {}

	~Heap() {
		delete[] m_innerHeap;
	}

	size_t size() const {
		return m_size;
	}

	bool empty() const {
		return m_size == 0;
	}

	const T & get(size_t index) const {
		assert(index < m_size);
		return m_innerHeap[index];
	}

	const T & minimum() const {
		assert(m_size > 0);
		return m_innerHeap[0];
	}

	const T & maximum() const {
		switch (m_size) {
		case 0: assert(false);
		case 1: return m_innerHeap[0];
		case 2: return m_innerHeap[1];
		default: /* >= 3*/ return m_innerHeap[2] < m_innerHeap[1] ? m_innerHeap[1] : m_innerHeap[2];
		}
	}

	void add(const T & value) {
		push_back(value);
		up(m_size - 1);
	}

	void clear() {
		m_size = 0;
	}

	void remove_min() {
		if (m_size == 0)
			return;
		--m_size;
		m_innerHeap[0] = m_innerHeap[m_size];
		down(0);
	}

	void remove_max() {
		if (m_size == 0)
			return;
		--m_size;
		if (m_size > 2) {
			if (m_innerHeap[2] < m_innerHeap[1]) {
				m_innerHeap[1] = m_innerHeap[m_size];
				down(1);
			}
			else {
				m_innerHeap[2] = m_innerHeap[m_size];
				down(2);
			}
		}
	}

	void remove(size_t index) {
		assert(index < m_size);
		--m_size;
		m_innerHeap[index] = m_innerHeap[m_size];
		down(index);
	}

	void updateValue(const T & newValue, size_t index) {
		assert(index < m_size);
		if (m_innerHeap[index] < newValue)
			up(index);
		else
			down(index);
	}
private:

	inline bool less(size_t index1, size_t index2, size_t level) {
		return (((level % 2) == 1 && m_innerHeap[index1] < m_innerHeap[index2]) ||
			((level % 2) == 0 && m_innerHeap[index2] < m_innerHeap[index1]));
	}

	static inline size_t log2(size_t x) {
		unsigned long index;
		_BitScanReverse64(&index, x); // always x > 0
		return index;
	}

	inline size_t getLevel(size_t index) {
		return log2(index + 1);
	}

	inline bool minLevel(size_t index) {
		return getLevel(index) % 2 == 0;
	}

	void down(size_t index) {
		if (minLevel(index))
			down_min(index);
		else
			down_max(index);
	}

	void down_min(size_t index) {
		do {
			size_t son = ((index + 1) << 1) - 1;
			if (son >= m_size)
				break;
			size_t grandchild = ((son + 1) << 1) - 1;
			size_t smallest = son;
			// Check the smallest among children
			++son;
			if (son < m_size && m_innerHeap[son] < m_innerHeap[smallest])
				smallest = son;
			// Check the smallest among grand-children
			for (size_t barrier = std::min(grandchild + 4, m_size); grandchild < barrier; ++grandchild) {
				if (m_innerHeap[grandchild] < m_innerHeap[smallest])
					smallest = grandchild;
			}
			if (son < smallest) { // Grandchild is the smallest
				if (m_innerHeap[smallest] < m_innerHeap[index]) {
					std::swap(m_innerHeap[smallest], m_innerHeap[index]);
					size_t parent = (smallest - 1) >> 1;
					if (m_innerHeap[parent] < m_innerHeap[smallest]) {
						std::swap(m_innerHeap[smallest], m_innerHeap[parent]);
					}
					index = smallest;
				}
				else
					break;
			}
			else { // Son is the smallest
				if (m_innerHeap[smallest] < m_innerHeap[index]) {
					std::swap(m_innerHeap[smallest], m_innerHeap[index]);
				}
				break;
			}
		} while (true);
	}

	void down_max(size_t index) {
		do {
			size_t son = ((index + 1) << 1) - 1;
			if (son >= m_size)
				break;
			size_t grandchild = ((son + 1) << 1) - 1;
			size_t largest = son;
			// Check the largest among children
			++son;
			if (son < m_size && m_innerHeap[largest] < m_innerHeap[son])
				largest = son;
			// Check the largest among grand-children
			for (size_t barrier = std::min(grandchild + 4, m_size); grandchild < barrier; ++grandchild) {
				if (m_innerHeap[largest] < m_innerHeap[grandchild])
					largest = grandchild;
			}
			if (son < largest) { // Grandchild is the largest
				if (m_innerHeap[index] < m_innerHeap[largest]) {
					std::swap(m_innerHeap[largest], m_innerHeap[index]);
					size_t parent = (largest - 1) >> 1;
					if (m_innerHeap[largest] < m_innerHeap[parent]) {
						std::swap(m_innerHeap[largest], m_innerHeap[parent]);
					}
					index = largest;
				}
				else
					break;
			}
			else { // Son is the largest
				if (m_innerHeap[index] < m_innerHeap[largest]) {
					std::swap(m_innerHeap[largest], m_innerHeap[index]);
				}
				break;
			}
		} while (true);
	}
	
	void up(size_t index) {
		size_t father = (index - 1) >> 1;
		if (minLevel(index)) {
			if (index > 0 && m_innerHeap[father] < m_innerHeap[index]) {
				std::swap(m_innerHeap[father], m_innerHeap[index]);
				up_max(father);
			}
			else
				up_min(index);
		}
		else {
			if (index > 0 && m_innerHeap[index] < m_innerHeap[father]) {
				std::swap(m_innerHeap[father], m_innerHeap[index]);
				up_min(father);
			}
			else
				up_max(index);
		}
	}

	void up_min(size_t index) {
		size_t grandparent = (((index - 1) >> 1) - 1) >> 1;
		while (index > 2 && m_innerHeap[index] < m_innerHeap[grandparent]) { // Has grand-parent
			std::swap(m_innerHeap[index], m_innerHeap[grandparent]);
			index = grandparent;
			grandparent = (((index - 1) >> 1) - 1) >> 1;
		}
	}

	void up_max(size_t index) {
		size_t grandparent = (((index - 1) >> 1) - 1) >> 1;
		while (index > 2 && m_innerHeap[grandparent] < m_innerHeap[index]) { // Has grand-parent
			std::swap(m_innerHeap[index], m_innerHeap[grandparent]);
			index = grandparent;
			grandparent = (((index - 1) >> 1) - 1) >> 1;
		}
	}

	void push_back(const T & value) {
		if (m_size == m_allocated) {
			size_t n = std::max((size_t)50, m_size * 2);
			T * temp = new T[n];
			memcpy(temp, m_innerHeap, m_size * sizeof(T));
			delete[] m_innerHeap;
			m_innerHeap = temp;
			m_allocated = n;
		}
		m_innerHeap[m_size++] = value;
	}
	T * m_innerHeap;
	size_t m_size, m_allocated;
};