#include <cstdio>
#include <fstream>
#include <string>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <map>
#include <set>

class GTFEntry {
public:
  std::string seqname;
  std::string source;
  std::string feature;
  int start;
  int end;
  double score;
  char strand;
  char frame;
  std::string attribute;

  int parse_gtf_line(const std::string& line)
  {
    std::size_t startpos = 0;
    std::size_t endpos = line.find('\t');
    seqname = line.substr(startpos, endpos-startpos);
    startpos = endpos + 1;

    endpos = line.find('\t', startpos);
    source = line.substr(startpos, endpos-startpos);
    startpos = endpos + 1;

    endpos = line.find('\t', startpos);
    feature = line.substr(startpos, endpos-startpos);
    startpos = endpos + 1;

    endpos = line.find('\t', startpos);
    start = std::atoi(line.substr(startpos, endpos-startpos).c_str());
    startpos = endpos + 1;

    endpos = line.find('\t', startpos);
    end = std::atoi(line.substr(startpos, endpos-startpos).c_str());
    startpos = endpos + 1;

    endpos = line.find('\t', startpos);
    std::string tmp = line.substr(startpos, endpos-startpos);
    if (tmp.compare(std::string(1, '.')) != 0)
      score = std::atoi(tmp.c_str());
    else
      score = 0.0;
    startpos = endpos + 1;

    endpos = line.find('\t', startpos);
    strand = line.at(startpos);
    startpos = endpos + 1;

    endpos = line.find('\t', startpos);
    frame = line.at(startpos);
    startpos = endpos + 1;

    endpos = line.find('\t', startpos);
    attribute = line.substr(startpos, endpos-startpos);
    startpos = endpos + 1;

    return 0;
  }

  void pretty_print(std::ostream& s)
  {
    s << "seqname: " << seqname << std::endl;
    s << "source: " << source << std::endl;
    s << "feature: " << feature << std::endl;
    s << "start: " << start << std::endl;
    s << "end: " << end << std::endl;
    s << "score: " << score << std::endl;
    s << "strand: " << strand << std::endl;
    s << "frame: " << frame << std::endl;
    s << "attribute: " << attribute << std::endl;
  }
};


void read_gtf(std::string file, std::vector<GTFEntry *>& gtf,
	      std::multimap<std::string, GTFEntry *>& gtf_per_chromosome,
	      std::set<std::string>& chromosomes)
{
  GTFEntry *entry;
  std::FILE *input;
  char buf[4096];
  std::string tmp;
  std::string line;

  tmp = "gunzip -c '";
  tmp += file;
  tmp += "'";
  input = popen(tmp.c_str(), "r");
  while (std::fgets(buf, 4096, input) != NULL) {
    line = std::string(buf);
    if (line.at(0) != '#') {
      entry = new GTFEntry();
      entry->parse_gtf_line(line);

      if (entry->feature.compare("gene") == 0) {
	gtf.push_back(entry);
	gtf_per_chromosome.insert(std::pair<std::string, GTFEntry *>(entry->seqname, entry));
	chromosomes.insert(entry->seqname);
      } else {
	delete entry;
      }
    }
  }
}



int main(int argc,char* argv[])
{
  std::vector<GTFEntry *> gtf = std::vector<GTFEntry *>();
  std::multimap<std::string, GTFEntry *> gtf_per_chromosome =
    std::multimap<std::string, GTFEntry *>();
  std::set<std::string> chromosomes = std::set<std::string>();

  read_gtf(argv[1], gtf, gtf_per_chromosome, chromosomes);

  //gtf[0]->pretty_print(std::cout);
  std::cout << gtf.size() << " genes found." << std::endl;
  for (std::set<std::string>::const_iterator i=chromosomes.begin();
       i != chromosomes.end(); i++) {
    std::cout << *i << ": " << gtf_per_chromosome.count(*i) 
	      << " genes." << std::endl;
  }

  return 0;
}
