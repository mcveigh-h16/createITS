perl ./check-cmsearch-v-cmscan.pl testfiles/1.cmsearch.tblout testfiles/1.cmscan.noclan.tblout
perl ./check-cmsearch-v-cmscan.pl testfiles/2.cmsearch.tblout testfiles/2.cmscan.noclan.tblout
perl ./check-cmsearch-v-cmscan.pl testfiles/3.cmsearch.tblout testfiles/3.cmscan.noclan.tblout

perl ./check-cmsearch-v-cmscan.pl --clanin ribo.claninfo testfiles/1.cmsearch.tblout testfiles/1.cmscan.clan.tblout
perl ./check-cmsearch-v-cmscan.pl --clanin ribo.claninfo testfiles/2.cmsearch.tblout testfiles/2.cmscan.clan.tblout
perl ./check-cmsearch-v-cmscan.pl --clanin ribo.claninfo testfiles/3.cmsearch.tblout testfiles/3.cmscan.clan.tblout
