// Prints what `istringstream i(line); i >> M >> k >> roff >> bc;` does for each input line
// (M, k, roff, bc are `unsigned int`, initialised to sentinels, the stream's default flags).
// Input: one line per row, with the escapes \t \n \r \v \f \0 \\ and \xHH. Output: tab separated
// `line  M  k  roff  bc  state`; a value is `-` when the stream left the sentinel in place;
// state is the stream's bits: g (good), e (eof), f (fail), b (bad).
#include <cstdio>
#include <iostream>
#include <sstream>
#include <string>

static std::string unescape(const std::string &t) {
    std::string out;
    for (size_t i = 0; i < t.size(); ++i) {
        if (t[i] != '\\') { out += t[i]; continue; }
        char c = t[++i];
        switch (c) {
        case 't': out += '\t'; break;
        case 'r': out += '\r'; break;
        case 'n': out += '\n'; break;
        case 'v': out += '\v'; break;
        case 'f': out += '\f'; break;
        case '0': out += '\0'; break;
        case '\\': out += '\\'; break;
        case 'x': out += (char)std::stoi(t.substr(i + 1, 2), nullptr, 16); i += 2; break;
        default: out += c;
        }
    }
    return out;
}

static std::string val(unsigned int v, unsigned int sentinel) {
    return v == sentinel ? "-" : std::to_string(v);
}

int main() {
    const unsigned int S = 0xA5A5A5A5u;
    std::string row;
    while (std::getline(std::cin, row)) {
        std::string line = unescape(row);
        unsigned int M = S, k = S, roff = S, bc = S;
        std::istringstream i(line);
        i >> M >> k >> roff >> bc;
        std::string st;
        if (i.good()) st += 'g';
        if (i.eof()) st += 'e';
        if (i.fail()) st += 'f';
        if (i.bad()) st += 'b';
        std::cout << row << '\t' << val(M, S) << '\t' << val(k, S) << '\t' << val(roff, S) << '\t'
                  << val(bc, S) << '\t' << st << '\n';
    }
    return 0;
}
