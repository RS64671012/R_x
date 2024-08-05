#pragma once
#include <iostream>
#include <vector>
#include <algorithm>
#include <limits>
#include <cmath>


struct Term {
	float coef;
	unsigned int exp;

	Term(float _coef = 0, unsigned int _exp = 0) {
		coef = _coef;
		exp = _exp;
	}

	Term operator- () const {
		return Term(-coef, exp);
	}

	Term operator* (const Term& other) const {
		return Term(coef * other.coef, exp + other.exp);
	}

	Term operator* (const float& other) const {
		return Term(coef * other, exp);
	}

	Term operator<< (const unsigned int& other) const {
		float c = coef;
		unsigned e = exp;
		if (c == 0 && e == 0) {
			c = float(std::numeric_limits<float>::quiet_NaN());
		}
		for (unsigned int i = 0; i < other; i++) {
			e += 1;
			c /= e;
		}
		return Term(c, e);
	}

	Term operator>> (const unsigned int& other) const {
		float c = coef;
		unsigned e = exp;
		for (unsigned int i = 0; i < other; i++) {
			if (e == 0) {
				c = 0;
				break;
			}
			c *= e;
			e -= 1;
		}
		return Term(c, e);
	}

	float at(const float& value) const {
		return coef * std::pow(value, exp);
	}

	float operator[] (const float& other) const {
		return at(other);
	}

	float evaluated(const float& start, const float& end) const {
		return at(end) - at(start);
	}

	void print(bool wrap = true) {
		float c = coef;
		if (std::isnan(c)) {
			std::cout << "C";
		}
		else if (c < 0) {
			c = -c;
			std::cout << "-";
		}
		else if (c != 1 || exp == 0) {
			std::cout << c;
		}
		if (exp == 1) {
			std::cout << "x";
		}
		else if (coef && exp) {
			std::cout << "x^" << exp;
		}
		if (wrap) {
			std::cout << std::endl;
		}
	}
};

static Term constant = Term(float(std::numeric_limits<float>::quiet_NaN()), 0);

class Polynomial {
private:
	std::vector<Term> terms;

	void combine() {
		std::sort(terms.begin(), terms.end(),
			[](const Term& a, const Term& b) {
				return a.exp > b.exp;
			});
		for (int i = 1; i < terms.size(); ) {
			if (terms[i].exp != terms[i - 1].exp) {
				i++;
				continue;
			}
			terms[i - 1].coef += terms[i].coef;
			terms.erase(terms.begin() + i);
		}
	}
	void zero_remove() {
		combine();
		for (int i = 0; i < terms.size() - 1; ) {
			if (terms[i].coef == 0) {
				terms.erase(terms.begin() + i);
				continue;
			}
			i++;
		}
	}

public:
	Polynomial(std::vector<Term> _terms = { Term() }) {
		terms = _terms;
		combine();
	}

	Polynomial(Term _term) {
		terms = { _term };
	}

	void fill() {
		for (int i = 0; i < terms.size(); i++) {
			addTerm(0, i);
		}
	}

	void addTerm(float coef, unsigned int exp) {
		terms.push_back(Term(coef, exp));
		combine();
	}

	void addTerm(Term term) {
		terms.push_back(term);
		combine();
	}

	Polynomial operator- () const {
		Polynomial result;
		for (int i = 0; i < terms.size(); i++) {
			result.addTerm(-terms[i]);
		}
		return result;
	}

	Polynomial operator+ (const Polynomial& other) const {
		Polynomial result;
		for (int i = 0; i < terms.size(); i++) {
			result.addTerm(terms[i]);
		}
		for (int j = 0; j < other.terms.size(); j++) {
			result.addTerm(other.terms[j]);
		}
		result.combine();
		return result;
	}

	Polynomial operator- (const Polynomial& other) const {
		Polynomial result;
		for (int i = 0; i < terms.size(); i++) {
			result.addTerm(terms[i]);
		}
		for (int j = 0; j < other.terms.size(); j++) {
			result.addTerm(-other.terms[j].coef, other.terms[j].exp);
		}
		result.combine();
		return result;
	}

	Polynomial operator* (const Polynomial& other) const {
		Polynomial result;
		for (int i = 0; i < terms.size(); i++) {
			for (int j = 0; j < other.terms.size(); j++) {
				result.addTerm(terms[i] * other.terms[j]);
			}
		}
		result.combine();
		return result;
	}

	Polynomial operator+ (const Term& other) const {
		Polynomial result = Polynomial(terms);
		result.addTerm(other);
		return result;
	}

	Polynomial operator- (const Term& other) const {
		Polynomial result;
		result.addTerm(-other.coef, other.exp);
		return result;
	}

	Polynomial operator* (const Term& other) const {
		Polynomial result;
		for (int i = 0; i < terms.size(); i++) {
			result.addTerm(terms[i] * other);
		}
		return result;
	}

	Polynomial operator/ (const Polynomial& other) const {
		Polynomial copy = Polynomial(terms), c_other = Polynomial(other), result;
		copy.zero_remove(); c_other.zero_remove();
		int exp_max = copy.terms[0].exp - c_other.terms[0].exp;
		copy.fill(); c_other.fill();
		if (exp_max < 0) {
			return Polynomial();
		}
		for (unsigned int shift = 0; shift <= exp_max; shift++) {
			float val = copy.terms[shift].coef / c_other.terms[0].coef;
			result.addTerm(val, exp_max - shift);
			for (int i = 0; i < c_other.terms.size(); i++) {
				copy.terms[i + shift].coef -= c_other.terms[i].coef * val;
			}
		}
		return result;
	}

	Polynomial operator% (const Polynomial& other) const {
		Polynomial copy = Polynomial(terms), c_other = Polynomial(other);
		copy.zero_remove(); c_other.zero_remove();
		int exp_max = copy.terms[0].exp - c_other.terms[0].exp;
		copy.fill(); c_other.fill();
		if (exp_max < 0) {
			return Polynomial(terms);
		}
		for (unsigned int shift = 0; shift <= exp_max; shift++) {
			float val = copy.terms[shift].coef / c_other.terms[0].coef;
			for (int i = 0; i < c_other.terms.size(); i++) {
				copy.terms[i + shift].coef -= c_other.terms[i].coef * val;
			}
		}
		return copy;
	}

	Polynomial operator<< (const unsigned int& other) const {
		Polynomial result = Polynomial(terms);
		for (int i = 0; i < terms.size(); i++) {
			result.terms[i] = result.terms[i] << other;
		}
		result.addTerm(constant);
		return result;
	}

	Polynomial operator>> (const unsigned int& other) const {
		Polynomial result = Polynomial(terms);
		for (int i = 0; i < terms.size(); i++) {
			result.terms[i] = result.terms[i] >> other;
		}
		return Polynomial(result);
	}

	Polynomial integral() const {
		return operator<< (1);
	}

	float integral(const float& start, const float& end) const {
		Polynomial result = Polynomial(terms);
		result.combine();
		if (result.terms[terms.size() - 1].coef == 0) {
			result.terms.pop_back();
		}
		for (int i = 0; i < result.terms.size(); i++) {
			result.terms[i] = result.terms[i] << 1;
		}
		return result.evaluated(start, end);
	}

	float at(const float& value) const {
		float result = 0;
		for (int i = 0; i < terms.size(); i++) {
			result += terms[i].at(value);
		}
		return result;
	}

	Polynomial operator& (const Polynomial& other) const {
		Polynomial copy = Polynomial(terms), c_other = Polynomial(other);
		while (c_other.terms[0].exp != 0) {
			copy = copy % c_other;
			std::vector<Polynomial> list = { c_other, copy };
			copy = list[0], c_other = list[1];
			copy.zero_remove(); c_other.zero_remove();
		}
		if (c_other.terms[0].coef != 0) {
			return Polynomial(Term(1));
		}
		return copy * Term(1 / other.terms[0].coef);
	}

	Polynomial gcd(const Polynomial& value) const {
		return operator&(value);
	}

	Polynomial operator| (const Polynomial& other) const {
		Polynomial copy = Polynomial(terms);
		return (other * copy) / (copy & other);
	}

	Polynomial lcm(const Polynomial& value) const {
		return operator|(value);
	}

	Polynomial operator^ (const Polynomial& other) const {
		Polynomial copy = Polynomial(terms);
		return (other | copy) / (copy & other);
	}

	float operator[] (const float& other) const {
		return at(other);
	}

	float evaluated(const float& start, const float& end) const {
		return at(end) - at(start);
	}

	void print() {
		bool with_sign = false;
		for (int i = 0; i < terms.size(); i++) {
			if (terms[i].coef == 0 && terms.size() != 1) { continue; }
			if (with_sign && (terms[i].coef > 0 || std::isnan(terms[i].coef))) { std::cout << "+"; }
			terms[i].print(false);
			with_sign = true;
		}
		std::cout << std::endl;
	}
};