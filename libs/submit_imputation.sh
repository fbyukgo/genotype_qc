TOKEN="eyJjdHkiOiJ0ZXh0XC9wbGFpbiIsImFsZyI6IkhTMjU2In0.eyJtYWlsIjoiYnV5dWtnb2xmdXJrYW5AZ21haWwuY29tIiwiYXBpX2hhc2giOiJiWU1kQ29nQWZOb0x0Q1ZBNks1UkJHU051dno4dFQiLCJleHBpcmUiOjE3Mzk0NTI0MDk3NzMsIm5hbWUiOiJGdXJrYW4gQsO8ecO8a2fDtmwiLCJhcGkiOnRydWUsInVzZXJuYW1lIjoiZnVya2FuMDE1ODAifQ.zH-Yi85nWm-U1jDrEpMB8uuy7p0X8_opOrPSVdxo4-0"


curl https://imputationserver.helmholtz-munich.de/api/v2/jobs/submit/imputationserver@2.0.0 \
  -H "X-Auth-Token: $TOKEN" \
  -F "files=@/ictstr01/groups/itg/teams/kim-hellmuth/users/furkan.bueyuekgoel/qtl/geno/data/vcf/PLINK_M00982_plate1-4-updated-chr1_refchecked.vcf.gz" \
  -F "refpanel=1000g-phase-3-v5" \
  -F "population=eur"

