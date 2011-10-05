// Generated by the protocol buffer compiler.  DO NOT EDIT!
// source: pp.proto

#ifndef PROTOBUF_pp_2eproto__INCLUDED
#define PROTOBUF_pp_2eproto__INCLUDED

#include <string>

#include <google/protobuf/stubs/common.h>

#if GOOGLE_PROTOBUF_VERSION < 2004000
#error This file was generated by a newer version of protoc which is
#error incompatible with your Protocol Buffer headers.  Please update
#error your headers.
#endif
#if 2004001 < GOOGLE_PROTOBUF_MIN_PROTOC_VERSION
#error This file was generated by an older version of protoc which is
#error incompatible with your Protocol Buffer headers.  Please
#error regenerate this file with a newer version of protoc.
#endif

#include <google/protobuf/generated_message_util.h>
#include <google/protobuf/repeated_field.h>
#include <google/protobuf/extension_set.h>
#include <google/protobuf/generated_message_reflection.h>
// @@protoc_insertion_point(includes)

// Internal implementation detail -- do not call these.
void  protobuf_AddDesc_pp_2eproto();
void protobuf_AssignDesc_pp_2eproto();
void protobuf_ShutdownFile_pp_2eproto();

class PolyPhen2Buffer;

enum PolyPhen2Buffer_pred_t {
  PolyPhen2Buffer_pred_t_UNKNOWN = 0,
  PolyPhen2Buffer_pred_t_BENIGN = 1,
  PolyPhen2Buffer_pred_t_POSS = 2,
  PolyPhen2Buffer_pred_t_PROB = 3
};
bool PolyPhen2Buffer_pred_t_IsValid(int value);
const PolyPhen2Buffer_pred_t PolyPhen2Buffer_pred_t_pred_t_MIN = PolyPhen2Buffer_pred_t_UNKNOWN;
const PolyPhen2Buffer_pred_t PolyPhen2Buffer_pred_t_pred_t_MAX = PolyPhen2Buffer_pred_t_PROB;
const int PolyPhen2Buffer_pred_t_pred_t_ARRAYSIZE = PolyPhen2Buffer_pred_t_pred_t_MAX + 1;

const ::google::protobuf::EnumDescriptor* PolyPhen2Buffer_pred_t_descriptor();
inline const ::std::string& PolyPhen2Buffer_pred_t_Name(PolyPhen2Buffer_pred_t value) {
  return ::google::protobuf::internal::NameOfEnum(
    PolyPhen2Buffer_pred_t_descriptor(), value);
}
inline bool PolyPhen2Buffer_pred_t_Parse(
    const ::std::string& name, PolyPhen2Buffer_pred_t* value) {
  return ::google::protobuf::internal::ParseNamedEnum<PolyPhen2Buffer_pred_t>(
    PolyPhen2Buffer_pred_t_descriptor(), name, value);
}
// ===================================================================

class PolyPhen2Buffer : public ::google::protobuf::Message {
 public:
  PolyPhen2Buffer();
  virtual ~PolyPhen2Buffer();
  
  PolyPhen2Buffer(const PolyPhen2Buffer& from);
  
  inline PolyPhen2Buffer& operator=(const PolyPhen2Buffer& from) {
    CopyFrom(from);
    return *this;
  }
  
  inline const ::google::protobuf::UnknownFieldSet& unknown_fields() const {
    return _unknown_fields_;
  }
  
  inline ::google::protobuf::UnknownFieldSet* mutable_unknown_fields() {
    return &_unknown_fields_;
  }
  
  static const ::google::protobuf::Descriptor* descriptor();
  static const PolyPhen2Buffer& default_instance();
  
  void Swap(PolyPhen2Buffer* other);
  
  // implements Message ----------------------------------------------
  
  PolyPhen2Buffer* New() const;
  void CopyFrom(const ::google::protobuf::Message& from);
  void MergeFrom(const ::google::protobuf::Message& from);
  void CopyFrom(const PolyPhen2Buffer& from);
  void MergeFrom(const PolyPhen2Buffer& from);
  void Clear();
  bool IsInitialized() const;
  
  int ByteSize() const;
  bool MergePartialFromCodedStream(
      ::google::protobuf::io::CodedInputStream* input);
  void SerializeWithCachedSizes(
      ::google::protobuf::io::CodedOutputStream* output) const;
  ::google::protobuf::uint8* SerializeWithCachedSizesToArray(::google::protobuf::uint8* output) const;
  int GetCachedSize() const { return _cached_size_; }
  private:
  void SharedCtor();
  void SharedDtor();
  void SetCachedSize(int size) const;
  public:
  
  ::google::protobuf::Metadata GetMetadata() const;
  
  // nested types ----------------------------------------------------
  
  typedef PolyPhen2Buffer_pred_t pred_t;
  static const pred_t UNKNOWN = PolyPhen2Buffer_pred_t_UNKNOWN;
  static const pred_t BENIGN = PolyPhen2Buffer_pred_t_BENIGN;
  static const pred_t POSS = PolyPhen2Buffer_pred_t_POSS;
  static const pred_t PROB = PolyPhen2Buffer_pred_t_PROB;
  static inline bool pred_t_IsValid(int value) {
    return PolyPhen2Buffer_pred_t_IsValid(value);
  }
  static const pred_t pred_t_MIN =
    PolyPhen2Buffer_pred_t_pred_t_MIN;
  static const pred_t pred_t_MAX =
    PolyPhen2Buffer_pred_t_pred_t_MAX;
  static const int pred_t_ARRAYSIZE =
    PolyPhen2Buffer_pred_t_pred_t_ARRAYSIZE;
  static inline const ::google::protobuf::EnumDescriptor*
  pred_t_descriptor() {
    return PolyPhen2Buffer_pred_t_descriptor();
  }
  static inline const ::std::string& pred_t_Name(pred_t value) {
    return PolyPhen2Buffer_pred_t_Name(value);
  }
  static inline bool pred_t_Parse(const ::std::string& name,
      pred_t* value) {
    return PolyPhen2Buffer_pred_t_Parse(name, value);
  }
  
  // accessors -------------------------------------------------------
  
  // required string transcript_name = 1;
  inline bool has_transcript_name() const;
  inline void clear_transcript_name();
  static const int kTranscriptNameFieldNumber = 1;
  inline const ::std::string& transcript_name() const;
  inline void set_transcript_name(const ::std::string& value);
  inline void set_transcript_name(const char* value);
  inline void set_transcript_name(const char* value, size_t size);
  inline ::std::string* mutable_transcript_name();
  inline ::std::string* release_transcript_name();
  
  // required string protein_name = 2;
  inline bool has_protein_name() const;
  inline void clear_protein_name();
  static const int kProteinNameFieldNumber = 2;
  inline const ::std::string& protein_name() const;
  inline void set_protein_name(const ::std::string& value);
  inline void set_protein_name(const char* value);
  inline void set_protein_name(const char* value, size_t size);
  inline ::std::string* mutable_protein_name();
  inline ::std::string* release_protein_name();
  
  // repeated int32 position = 3 [packed = true];
  inline int position_size() const;
  inline void clear_position();
  static const int kPositionFieldNumber = 3;
  inline ::google::protobuf::int32 position(int index) const;
  inline void set_position(int index, ::google::protobuf::int32 value);
  inline void add_position(::google::protobuf::int32 value);
  inline const ::google::protobuf::RepeatedField< ::google::protobuf::int32 >&
      position() const;
  inline ::google::protobuf::RepeatedField< ::google::protobuf::int32 >*
      mutable_position();
  
  // repeated string reference = 4;
  inline int reference_size() const;
  inline void clear_reference();
  static const int kReferenceFieldNumber = 4;
  inline const ::std::string& reference(int index) const;
  inline ::std::string* mutable_reference(int index);
  inline void set_reference(int index, const ::std::string& value);
  inline void set_reference(int index, const char* value);
  inline void set_reference(int index, const char* value, size_t size);
  inline ::std::string* add_reference();
  inline void add_reference(const ::std::string& value);
  inline void add_reference(const char* value);
  inline void add_reference(const char* value, size_t size);
  inline const ::google::protobuf::RepeatedPtrField< ::std::string>& reference() const;
  inline ::google::protobuf::RepeatedPtrField< ::std::string>* mutable_reference();
  
  // repeated string alternate = 5;
  inline int alternate_size() const;
  inline void clear_alternate();
  static const int kAlternateFieldNumber = 5;
  inline const ::std::string& alternate(int index) const;
  inline ::std::string* mutable_alternate(int index);
  inline void set_alternate(int index, const ::std::string& value);
  inline void set_alternate(int index, const char* value);
  inline void set_alternate(int index, const char* value, size_t size);
  inline ::std::string* add_alternate();
  inline void add_alternate(const ::std::string& value);
  inline void add_alternate(const char* value);
  inline void add_alternate(const char* value, size_t size);
  inline const ::google::protobuf::RepeatedPtrField< ::std::string>& alternate() const;
  inline ::google::protobuf::RepeatedPtrField< ::std::string>* mutable_alternate();
  
  // repeated double score = 6 [packed = true];
  inline int score_size() const;
  inline void clear_score();
  static const int kScoreFieldNumber = 6;
  inline double score(int index) const;
  inline void set_score(int index, double value);
  inline void add_score(double value);
  inline const ::google::protobuf::RepeatedField< double >&
      score() const;
  inline ::google::protobuf::RepeatedField< double >*
      mutable_score();
  
  // repeated .PolyPhen2Buffer.pred_t prediction = 7 [packed = true];
  inline int prediction_size() const;
  inline void clear_prediction();
  static const int kPredictionFieldNumber = 7;
  inline ::PolyPhen2Buffer_pred_t prediction(int index) const;
  inline void set_prediction(int index, ::PolyPhen2Buffer_pred_t value);
  inline void add_prediction(::PolyPhen2Buffer_pred_t value);
  inline const ::google::protobuf::RepeatedField<int>& prediction() const;
  inline ::google::protobuf::RepeatedField<int>* mutable_prediction();
  
  // @@protoc_insertion_point(class_scope:PolyPhen2Buffer)
 private:
  inline void set_has_transcript_name();
  inline void clear_has_transcript_name();
  inline void set_has_protein_name();
  inline void clear_has_protein_name();
  
  ::google::protobuf::UnknownFieldSet _unknown_fields_;
  
  ::std::string* transcript_name_;
  ::std::string* protein_name_;
  ::google::protobuf::RepeatedField< ::google::protobuf::int32 > position_;
  mutable int _position_cached_byte_size_;
  ::google::protobuf::RepeatedPtrField< ::std::string> reference_;
  ::google::protobuf::RepeatedPtrField< ::std::string> alternate_;
  ::google::protobuf::RepeatedField< double > score_;
  mutable int _score_cached_byte_size_;
  ::google::protobuf::RepeatedField<int> prediction_;
  mutable int _prediction_cached_byte_size_;
  
  mutable int _cached_size_;
  ::google::protobuf::uint32 _has_bits_[(7 + 31) / 32];
  
  friend void  protobuf_AddDesc_pp_2eproto();
  friend void protobuf_AssignDesc_pp_2eproto();
  friend void protobuf_ShutdownFile_pp_2eproto();
  
  void InitAsDefaultInstance();
  static PolyPhen2Buffer* default_instance_;
};
// ===================================================================


// ===================================================================

// PolyPhen2Buffer

// required string transcript_name = 1;
inline bool PolyPhen2Buffer::has_transcript_name() const {
  return (_has_bits_[0] & 0x00000001u) != 0;
}
inline void PolyPhen2Buffer::set_has_transcript_name() {
  _has_bits_[0] |= 0x00000001u;
}
inline void PolyPhen2Buffer::clear_has_transcript_name() {
  _has_bits_[0] &= ~0x00000001u;
}
inline void PolyPhen2Buffer::clear_transcript_name() {
  if (transcript_name_ != &::google::protobuf::internal::kEmptyString) {
    transcript_name_->clear();
  }
  clear_has_transcript_name();
}
inline const ::std::string& PolyPhen2Buffer::transcript_name() const {
  return *transcript_name_;
}
inline void PolyPhen2Buffer::set_transcript_name(const ::std::string& value) {
  set_has_transcript_name();
  if (transcript_name_ == &::google::protobuf::internal::kEmptyString) {
    transcript_name_ = new ::std::string;
  }
  transcript_name_->assign(value);
}
inline void PolyPhen2Buffer::set_transcript_name(const char* value) {
  set_has_transcript_name();
  if (transcript_name_ == &::google::protobuf::internal::kEmptyString) {
    transcript_name_ = new ::std::string;
  }
  transcript_name_->assign(value);
}
inline void PolyPhen2Buffer::set_transcript_name(const char* value, size_t size) {
  set_has_transcript_name();
  if (transcript_name_ == &::google::protobuf::internal::kEmptyString) {
    transcript_name_ = new ::std::string;
  }
  transcript_name_->assign(reinterpret_cast<const char*>(value), size);
}
inline ::std::string* PolyPhen2Buffer::mutable_transcript_name() {
  set_has_transcript_name();
  if (transcript_name_ == &::google::protobuf::internal::kEmptyString) {
    transcript_name_ = new ::std::string;
  }
  return transcript_name_;
}
inline ::std::string* PolyPhen2Buffer::release_transcript_name() {
  clear_has_transcript_name();
  if (transcript_name_ == &::google::protobuf::internal::kEmptyString) {
    return NULL;
  } else {
    ::std::string* temp = transcript_name_;
    transcript_name_ = const_cast< ::std::string*>(&::google::protobuf::internal::kEmptyString);
    return temp;
  }
}

// required string protein_name = 2;
inline bool PolyPhen2Buffer::has_protein_name() const {
  return (_has_bits_[0] & 0x00000002u) != 0;
}
inline void PolyPhen2Buffer::set_has_protein_name() {
  _has_bits_[0] |= 0x00000002u;
}
inline void PolyPhen2Buffer::clear_has_protein_name() {
  _has_bits_[0] &= ~0x00000002u;
}
inline void PolyPhen2Buffer::clear_protein_name() {
  if (protein_name_ != &::google::protobuf::internal::kEmptyString) {
    protein_name_->clear();
  }
  clear_has_protein_name();
}
inline const ::std::string& PolyPhen2Buffer::protein_name() const {
  return *protein_name_;
}
inline void PolyPhen2Buffer::set_protein_name(const ::std::string& value) {
  set_has_protein_name();
  if (protein_name_ == &::google::protobuf::internal::kEmptyString) {
    protein_name_ = new ::std::string;
  }
  protein_name_->assign(value);
}
inline void PolyPhen2Buffer::set_protein_name(const char* value) {
  set_has_protein_name();
  if (protein_name_ == &::google::protobuf::internal::kEmptyString) {
    protein_name_ = new ::std::string;
  }
  protein_name_->assign(value);
}
inline void PolyPhen2Buffer::set_protein_name(const char* value, size_t size) {
  set_has_protein_name();
  if (protein_name_ == &::google::protobuf::internal::kEmptyString) {
    protein_name_ = new ::std::string;
  }
  protein_name_->assign(reinterpret_cast<const char*>(value), size);
}
inline ::std::string* PolyPhen2Buffer::mutable_protein_name() {
  set_has_protein_name();
  if (protein_name_ == &::google::protobuf::internal::kEmptyString) {
    protein_name_ = new ::std::string;
  }
  return protein_name_;
}
inline ::std::string* PolyPhen2Buffer::release_protein_name() {
  clear_has_protein_name();
  if (protein_name_ == &::google::protobuf::internal::kEmptyString) {
    return NULL;
  } else {
    ::std::string* temp = protein_name_;
    protein_name_ = const_cast< ::std::string*>(&::google::protobuf::internal::kEmptyString);
    return temp;
  }
}

// repeated int32 position = 3 [packed = true];
inline int PolyPhen2Buffer::position_size() const {
  return position_.size();
}
inline void PolyPhen2Buffer::clear_position() {
  position_.Clear();
}
inline ::google::protobuf::int32 PolyPhen2Buffer::position(int index) const {
  return position_.Get(index);
}
inline void PolyPhen2Buffer::set_position(int index, ::google::protobuf::int32 value) {
  position_.Set(index, value);
}
inline void PolyPhen2Buffer::add_position(::google::protobuf::int32 value) {
  position_.Add(value);
}
inline const ::google::protobuf::RepeatedField< ::google::protobuf::int32 >&
PolyPhen2Buffer::position() const {
  return position_;
}
inline ::google::protobuf::RepeatedField< ::google::protobuf::int32 >*
PolyPhen2Buffer::mutable_position() {
  return &position_;
}

// repeated string reference = 4;
inline int PolyPhen2Buffer::reference_size() const {
  return reference_.size();
}
inline void PolyPhen2Buffer::clear_reference() {
  reference_.Clear();
}
inline const ::std::string& PolyPhen2Buffer::reference(int index) const {
  return reference_.Get(index);
}
inline ::std::string* PolyPhen2Buffer::mutable_reference(int index) {
  return reference_.Mutable(index);
}
inline void PolyPhen2Buffer::set_reference(int index, const ::std::string& value) {
  reference_.Mutable(index)->assign(value);
}
inline void PolyPhen2Buffer::set_reference(int index, const char* value) {
  reference_.Mutable(index)->assign(value);
}
inline void PolyPhen2Buffer::set_reference(int index, const char* value, size_t size) {
  reference_.Mutable(index)->assign(
    reinterpret_cast<const char*>(value), size);
}
inline ::std::string* PolyPhen2Buffer::add_reference() {
  return reference_.Add();
}
inline void PolyPhen2Buffer::add_reference(const ::std::string& value) {
  reference_.Add()->assign(value);
}
inline void PolyPhen2Buffer::add_reference(const char* value) {
  reference_.Add()->assign(value);
}
inline void PolyPhen2Buffer::add_reference(const char* value, size_t size) {
  reference_.Add()->assign(reinterpret_cast<const char*>(value), size);
}
inline const ::google::protobuf::RepeatedPtrField< ::std::string>&
PolyPhen2Buffer::reference() const {
  return reference_;
}
inline ::google::protobuf::RepeatedPtrField< ::std::string>*
PolyPhen2Buffer::mutable_reference() {
  return &reference_;
}

// repeated string alternate = 5;
inline int PolyPhen2Buffer::alternate_size() const {
  return alternate_.size();
}
inline void PolyPhen2Buffer::clear_alternate() {
  alternate_.Clear();
}
inline const ::std::string& PolyPhen2Buffer::alternate(int index) const {
  return alternate_.Get(index);
}
inline ::std::string* PolyPhen2Buffer::mutable_alternate(int index) {
  return alternate_.Mutable(index);
}
inline void PolyPhen2Buffer::set_alternate(int index, const ::std::string& value) {
  alternate_.Mutable(index)->assign(value);
}
inline void PolyPhen2Buffer::set_alternate(int index, const char* value) {
  alternate_.Mutable(index)->assign(value);
}
inline void PolyPhen2Buffer::set_alternate(int index, const char* value, size_t size) {
  alternate_.Mutable(index)->assign(
    reinterpret_cast<const char*>(value), size);
}
inline ::std::string* PolyPhen2Buffer::add_alternate() {
  return alternate_.Add();
}
inline void PolyPhen2Buffer::add_alternate(const ::std::string& value) {
  alternate_.Add()->assign(value);
}
inline void PolyPhen2Buffer::add_alternate(const char* value) {
  alternate_.Add()->assign(value);
}
inline void PolyPhen2Buffer::add_alternate(const char* value, size_t size) {
  alternate_.Add()->assign(reinterpret_cast<const char*>(value), size);
}
inline const ::google::protobuf::RepeatedPtrField< ::std::string>&
PolyPhen2Buffer::alternate() const {
  return alternate_;
}
inline ::google::protobuf::RepeatedPtrField< ::std::string>*
PolyPhen2Buffer::mutable_alternate() {
  return &alternate_;
}

// repeated double score = 6 [packed = true];
inline int PolyPhen2Buffer::score_size() const {
  return score_.size();
}
inline void PolyPhen2Buffer::clear_score() {
  score_.Clear();
}
inline double PolyPhen2Buffer::score(int index) const {
  return score_.Get(index);
}
inline void PolyPhen2Buffer::set_score(int index, double value) {
  score_.Set(index, value);
}
inline void PolyPhen2Buffer::add_score(double value) {
  score_.Add(value);
}
inline const ::google::protobuf::RepeatedField< double >&
PolyPhen2Buffer::score() const {
  return score_;
}
inline ::google::protobuf::RepeatedField< double >*
PolyPhen2Buffer::mutable_score() {
  return &score_;
}

// repeated .PolyPhen2Buffer.pred_t prediction = 7 [packed = true];
inline int PolyPhen2Buffer::prediction_size() const {
  return prediction_.size();
}
inline void PolyPhen2Buffer::clear_prediction() {
  prediction_.Clear();
}
inline ::PolyPhen2Buffer_pred_t PolyPhen2Buffer::prediction(int index) const {
  return static_cast< ::PolyPhen2Buffer_pred_t >(prediction_.Get(index));
}
inline void PolyPhen2Buffer::set_prediction(int index, ::PolyPhen2Buffer_pred_t value) {
  GOOGLE_DCHECK(::PolyPhen2Buffer_pred_t_IsValid(value));
  prediction_.Set(index, value);
}
inline void PolyPhen2Buffer::add_prediction(::PolyPhen2Buffer_pred_t value) {
  GOOGLE_DCHECK(::PolyPhen2Buffer_pred_t_IsValid(value));
  prediction_.Add(value);
}
inline const ::google::protobuf::RepeatedField<int>&
PolyPhen2Buffer::prediction() const {
  return prediction_;
}
inline ::google::protobuf::RepeatedField<int>*
PolyPhen2Buffer::mutable_prediction() {
  return &prediction_;
}


// @@protoc_insertion_point(namespace_scope)

#ifndef SWIG
namespace google {
namespace protobuf {

template <>
inline const EnumDescriptor* GetEnumDescriptor< ::PolyPhen2Buffer_pred_t>() {
  return ::PolyPhen2Buffer_pred_t_descriptor();
}

}  // namespace google
}  // namespace protobuf
#endif  // SWIG

// @@protoc_insertion_point(global_scope)

#endif  // PROTOBUF_pp_2eproto__INCLUDED
