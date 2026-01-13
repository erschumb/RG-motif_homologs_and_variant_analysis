########################### PUT HERE API ACCESS AND POSSIBLE FILTERING?

# This example uses the GQL GraphQL client library.
#
# To install: pip3 install gql
#
# GQL is one popular Python GraphQL client, but there are others.
# See https://graphql.org/community/tools-and-libraries/?tags=python_client

from gql import gql, Client
from gql.transport.aiohttp import AIOHTTPTransport

transport = AIOHTTPTransport(url="https://gnomad.broadinstitute.org/api")
client = Client(transport=transport, fetch_schema_from_transport=True)

# For brevity, and to keep the focus on the Python code, we don't include every
# field from the raw query here.

query = gql(
    """
    query VariantsInGene {
      gene(gene_symbol: "MAP7D1", reference_genome: GRCh38) {
        variants(dataset: gnomad_r4) {
            variant_id
            pos
            rsids
            transcript_id
            transcript_version
            hgvs
            hgvsc
            hgvsp
            flags
            consequence
            exome {
                ac
                ac_hemi
                ac_hom
                an
                af
                populations {
                    id
                    ac
                    an
                    ac_hemi
                    ac_hom
                }
            }
        }
      }
    }
"""
)

# Execute the query on the transport


##### result = await client.execute_async(query)


# print(list(result.keys()))
# print(list(result["gene"].keys()))
for el in result["gene"]["variants"]:
    if el["pos"] == 36176314:
        print(el)
# print(result["gene"]["variants"])

#### WRITE THIS TO MAKE THE SAME OUTPUT AS THE ANALYSIS THROUGH THE BED FILES AND MERTS SCRIPT

