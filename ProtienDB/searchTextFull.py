from rcsbapi.search import TextQuery

# Search for structures associated with the phrase "Hemoglobin"
query = TextQuery(value="Aspirin")

# Execute the query by running it as a function
results = query()

# Results are returned as an iterator of result identifiers.
for rid in results:
    print(rid)

