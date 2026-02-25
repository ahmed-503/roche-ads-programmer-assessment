# question_4_genai/clinical_data_agent.py

from __future__ import annotations

import json
import os
import re
from dataclasses import dataclass
from typing import Any, Dict, List, Optional

import pandas as pd

# Optional LangChain/OpenAI imports (only used if available and API key is set)
try:
    from langchain_core.prompts import ChatPromptTemplate
    from langchain_core.output_parsers import StrOutputParser
    from langchain_openai import ChatOpenAI
    LANGCHAIN_AVAILABLE = True
except Exception:
    LANGCHAIN_AVAILABLE = False


@dataclass
class ParsedQuery:
    target_column: str
    filter_value: str


class ClinicalTrialDataAgent:
    """
    GenAI Clinical Data Assistant for AE dataset.
    
    Responsibilities:
    1) Understand schema (column meanings).
    2) Convert natural language questions into structured JSON:
       - target_column
       - filter_value
    3) Execute pandas filter on AE dataframe.
    4) Return unique subject count + list of subject IDs.
    """

    def __init__(
        self,
        ae_df: pd.DataFrame,
        use_llm: bool = True,
        model_name: str = "gpt-4o-mini",
        openai_api_key: Optional[str] = None,
    ) -> None:
        self.ae_df = ae_df.copy()
        self.use_llm = use_llm
        self.model_name = model_name
        self.openai_api_key = openai_api_key or os.getenv("OPENAI_API_KEY")

        # Minimal schema description required by prompt
        self.schema_description = """
Dataset: AE (Adverse Events)
Relevant columns:
- USUBJID: Unique Subject Identifier (use this to count unique subjects and return matching IDs)
- AETERM: Reported term for adverse event / condition (e.g., HEADACHE, NAUSEA)
- AESOC: System organ class / body system (e.g., CARDIAC DISORDERS, SKIN AND SUBCUTANEOUS TISSUE DISORDERS)
- AESEV: Severity / intensity of adverse event (e.g., MILD, MODERATE, SEVERE)

Task:
Given a user question, identify the most appropriate target column and extract the filter value.

Output JSON only:
{
  "target_column": "<AESEV|AETERM|AESOC>",
  "filter_value": "<string>"
}
""".strip()

        # Basic validation
        required_cols = {"USUBJID", "AETERM", "AESOC", "AESEV"}
        missing = required_cols - set(self.ae_df.columns)
        if missing:
            raise ValueError(f"AE dataset is missing required columns: {sorted(missing)}")

    # ----------------------------
    # Public interface
    # ----------------------------
    def ask(self, question: str) -> Dict[str, Any]:
        """
        Full flow:
        Prompt -> Parse -> Execute
        """
        parsed = self.parse_question(question)
        result = self.execute_filter(parsed)
        return {
            "question": question,
            "parsed_query": {
                "target_column": parsed.target_column,
                "filter_value": parsed.filter_value,
            },
            "result": result,
        }

    def parse_question(self, question: str) -> ParsedQuery:
        """
        Use LLM if available; otherwise fallback to mock/rule-based parser.
        """
        if self.use_llm and LANGCHAIN_AVAILABLE and self.openai_api_key:
            try:
                return self._parse_with_llm(question)
            except Exception as e:
                # Fall back gracefully if the LLM call fails
                print(f"[WARN] LLM parsing failed, falling back to mock parser. Error: {e}")

        return self._mock_parse(question)

    def execute_filter(self, parsed_query: ParsedQuery) -> Dict[str, Any]:
        """
        Apply pandas filtering and return unique USUBJIDs + count.
        """
        col = parsed_query.target_column
        val = str(parsed_query.filter_value).strip()

        if col not in self.ae_df.columns:
            raise ValueError(f"Invalid target_column '{col}'. Available columns include: {list(self.ae_df.columns)}")

        # Case-insensitive exact/contains hybrid matching:
        # - Exact match first
        # - If no exact match, fallback to contains
        series = self.ae_df[col].astype(str)

        exact_mask = series.str.upper() == val.upper()
        if exact_mask.any():
            matched = self.ae_df[exact_mask].copy()
        else:
            contains_mask = series.str.contains(re.escape(val), case=False, na=False)
            matched = self.ae_df[contains_mask].copy()

        # Unique subjects
        subject_ids = (
            matched["USUBJID"]
            .dropna()
            .astype(str)
            .drop_duplicates()
            .sort_values()
            .tolist()
        )

        return {
            "target_column": col,
            "filter_value": val,
            "n_matching_rows": int(len(matched)),
            "n_unique_subjects": int(len(subject_ids)),
            "matching_usubjid": subject_ids,
        }

    # ----------------------------
    # LLM parser (LangChain + OpenAI)
    # ----------------------------
    def _parse_with_llm(self, question: str) -> ParsedQuery:
        """
        Uses LangChain + OpenAI to return structured JSON text, then parses it.
        """
        if not LANGCHAIN_AVAILABLE:
            raise RuntimeError("LangChain/OpenAI dependencies are not installed.")

        llm = ChatOpenAI(
            model=self.model_name,
            temperature=0,
            api_key=self.openai_api_key,
        )

        prompt = ChatPromptTemplate.from_messages(
            [
                (
                    "system",
                    "You are a clinical data query parser. "
                    "Map user intent to the correct AE dataset column and extract the filter value.\n\n"
                    f"{self.schema_description}\n\n"
                    "Rules:\n"
                    "- If user asks about severity/intensity -> AESEV\n"
                    "- If user asks about a condition/event term (e.g., headache, nausea) -> AETERM\n"
                    "- If user asks about a body system (e.g., cardiac, skin) -> AESOC\n"
                    "- Return valid JSON only, no markdown."
                ),
                ("human", "Question: {question}"),
            ]
        )

        chain = prompt | llm | StrOutputParser()
        raw = chain.invoke({"question": question})

        parsed = self._safe_json_parse(raw)
        target_column = parsed.get("target_column", "").strip().upper()
        filter_value = str(parsed.get("filter_value", "")).strip()

        if target_column not in {"AESEV", "AETERM", "AESOC"}:
            raise ValueError(f"LLM returned unsupported target_column: {target_column}")
        if not filter_value:
            raise ValueError("LLM returned empty filter_value")

        return ParsedQuery(target_column=target_column, filter_value=filter_value)

    # ----------------------------
    # Mock parser (works without API key)
    # ----------------------------
    def _mock_parse(self, question: str) -> ParsedQuery:
        """
        Rule-based/mock parser to satisfy the assessment if no API key is available.
        Still preserves Prompt -> Parse -> Execute logic.
        """
        q = question.strip()
        q_lower = q.lower()

        # 1) Severity / intensity -> AESEV
        if any(word in q_lower for word in ["severity", "intensity", "severe", "moderate", "mild"]):
            # Try to extract known severity values
            for sev in ["MILD", "MODERATE", "SEVERE"]:
                if sev.lower() in q_lower:
                    return ParsedQuery(target_column="AESEV", filter_value=sev)

            # If "severity" mentioned but value not obvious, capture a phrase after "of"
            # Example: "AEs of moderate severity"
            match = re.search(r"\b(?:of|with)\s+([a-zA-Z]+)\s+severity\b", q_lower)
            if match:
                return ParsedQuery(target_column="AESEV", filter_value=match.group(1).upper())

            # Fallback
            return ParsedQuery(target_column="AESEV", filter_value="MODERATE")

        # 2) Body system -> AESOC
        if any(word in q_lower for word in ["body system", "system organ class", "soc", "cardiac", "skin", "gastro", "nervous"]):
            # Common shortcuts
            if "cardiac" in q_lower:
                return ParsedQuery(target_column="AESOC", filter_value="CARDIAC")
            if "skin" in q_lower:
                return ParsedQuery(target_column="AESOC", filter_value="SKIN")
            if "gastro" in q_lower:
                return ParsedQuery(target_column="AESOC", filter_value="GASTRO")
            if "nervous" in q_lower:
                return ParsedQuery(target_column="AESOC", filter_value="NERVOUS")

            # Try quoted phrase
            quoted = self._extract_quoted_text(q)
            if quoted:
                return ParsedQuery(target_column="AESOC", filter_value=quoted)

            return ParsedQuery(target_column="AESOC", filter_value="CARDIAC")

        # 3) Condition/event term -> AETERM (default for specific events)
        # Try quoted phrase first
        quoted = self._extract_quoted_text(q)
        if quoted:
            return ParsedQuery(target_column="AETERM", filter_value=quoted)

        # Heuristic: look for "had X", "with X", "event X"
        event_match = re.search(
            r"(?:had|with|for|event|events of)\s+([A-Za-z][A-Za-z\s\-\/]+?)(?:\?|$)",
            q,
            flags=re.IGNORECASE,
        )
        if event_match:
            candidate = event_match.group(1).strip()
            candidate = re.sub(r"\b(adverse events?|aes?)\b", "", candidate, flags=re.IGNORECASE).strip()
            candidate = re.sub(r"\s+", " ", candidate)
            if candidate:
                return ParsedQuery(target_column="AETERM", filter_value=candidate)

        # Final fallback
        return ParsedQuery(target_column="AETERM", filter_value="HEADACHE")

    # ----------------------------
    # Helpers
    # ----------------------------
    @staticmethod
    def _extract_quoted_text(text: str) -> Optional[str]:
        match = re.search(r'["“](.+?)["”]', text)
        return match.group(1).strip() if match else None

    @staticmethod
    def _safe_json_parse(raw_text: str) -> Dict[str, Any]:
        """
        Robust JSON extraction in case the LLM returns extra text.
        """
        raw_text = raw_text.strip()

        # Try direct parse first
        try:
            return json.loads(raw_text)
        except json.JSONDecodeError:
            pass

        # Try to extract JSON object substring
        start = raw_text.find("{")
        end = raw_text.rfind("}")
        if start != -1 and end != -1 and end > start:
            return json.loads(raw_text[start:end + 1])

        raise ValueError(f"Could not parse JSON from LLM output: {raw_text}")


# ----------------------------------
# Utility: load AE data
# ----------------------------------
def load_ae_csv(csv_path: str) -> pd.DataFrame:
    """
    Load adae.csv (or ae-like file) into pandas.
    """
    df = pd.read_csv(csv_path)

    # Normalize column names in case lowercase file is provided
    df.columns = [c.upper() for c in df.columns]

    return df


# ----------------------------------
# Test script (3 example queries)
# ----------------------------------
if __name__ == "__main__":
    # Update this path to your actual file location
    # The assessment says "adae.csv (pharmaversesdtm::ae)" — use the CSV you export/provide.
    CSV_PATH = "adae.csv"

    ae = load_ae_csv(CSV_PATH)

    # Set use_llm=True if you have OPENAI_API_KEY and langchain packages installed.
    # Otherwise leave False to use the mock parser.
    agent = ClinicalTrialDataAgent(
        ae_df=ae,
        use_llm=False,  # change to True if using OpenAI/LangChain
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