// This tells Catch to provide a main() - only do this in one cpp file.
#define CATCH_CONFIG_MAIN

#include "catch.hpp"
#include "partition.hpp"



namespace pt {


TEST_CASE("Partition", "[partition]") {

 auto p = std::unique_ptr<Partition>(
   new Partition("test-data/five/RAxML_bestTree.five", "test-data/five/five.fasta"));

 std::cout << p->ToNewick() << std::endl;
 REQUIRE(p->branch_count() == 7);
 REQUIRE(-4159.020246 - p->FullTraversalLogLikelihood() < 1e-6);

}

}
