# question_4_genai/run_examples.py

import json
from clinical_data_agent import ClinicalTrialDataAgent, load_ae_csv

CSV_PATH = "adae.csv"

ae = load_ae_csv(CSV_PATH)

agent = ClinicalTrialDataAgent(
    ae_df=ae,
    use_llm=False,  # mock mode for assessment/demo
    model_name="gpt-4o-mini",
)

example_questions = [
    "Give me the subjects who had Adverse events of Moderate severity.",
    'Which subjects had the condition "Headache"?',
    "Show me subjects with cardiac adverse events.",
]

for i, q in enumerate(example_questions, start=1):
    print("\n" + "=" * 80)
    print(f"Example Query {i}: {q}")
    response = agent.ask(q)
    print(json.dumps(response, indent=2))